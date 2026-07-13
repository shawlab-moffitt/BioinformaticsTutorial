
#!/usr/bin/env Rscript

# Download and format GSE36133 for DRPPM-EASY / Easy-App
#
# Outputs:
#   GSE36133_expression.tsv
#   GSE36133_metadata.tsv
#   GSE36133_expression_probe_level.tsv
#   GSE36133_GEO_objects.rds
#
# Expression output format:
#   Gene    GSM886835    GSM886836 ...
#
# Metadata output format:
#   geo_accession    title    ...metadata columns...

# -------------------------------------------------------------------------
# Install required packages if necessary
# -------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c("GEOquery", "Biobase")

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

cran_packages <- c("dplyr", "tidyr", "readr", "tibble", "stringr")

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(stringr)
})

# -------------------------------------------------------------------------
# Settings
# -------------------------------------------------------------------------

gse_id <- "GSE36133"
output_dir <- "GSE36133_EasyApp"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Downloading ", gse_id, " from GEO...")

# This follows the Easy-App workflow:
# gset <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = TRUE)
gset_list <- GEOquery::getGEO(
  GEO = gse_id,
  GSEMatrix = TRUE,
  getGPL = TRUE,
  destdir = output_dir
)

if (length(gset_list) == 0) {
  stop("No ExpressionSet was returned for ", gse_id)
}

message(
  "Number of platform-specific ExpressionSets returned: ",
  length(gset_list)
)

# GSE36133 should contain GPL15308. Select it explicitly when available.
gset_names <- names(gset_list)
platform_index <- grep("GPL15308", gset_names)

if (length(platform_index) == 0) {
  warning(
    "GPL15308 was not found in the ExpressionSet names. ",
    "Using the first ExpressionSet."
  )
  platform_index <- 1
} else {
  platform_index <- platform_index[1]
}

gset <- gset_list[[platform_index]]

platform_id <- annotation(gset)

if (is.null(platform_id) || platform_id == "") {
  platform_id <- sub("^.*-", "", gset_names[platform_index])
}

message("Selected platform: ", platform_id)
message("Samples: ", ncol(exprs(gset)))
message("Features: ", nrow(exprs(gset)))

# -------------------------------------------------------------------------
# Extract the original expression matrix
# -------------------------------------------------------------------------

expression_matrix <- Biobase::exprs(gset)

expression_probe <- as.data.frame(
  expression_matrix,
  check.names = FALSE
) %>%
  rownames_to_column("ID")

write_tsv(
  expression_probe,
  file.path(output_dir, paste0(gse_id, "_expression_probe_level.tsv"))
)

# -------------------------------------------------------------------------
# Obtain platform annotation
# -------------------------------------------------------------------------

# Expression matrix
expression_matrix <- Biobase::exprs(gset)

# GPL15308 uses Brainarray HGU133Plus2_Hs_ENTREZG identifiers.
# Row names may look like "7157_at" or simply "7157".
entrez_ids <- rownames(expression_matrix)
entrez_ids <- sub("_at$", "", entrez_ids)

# Install annotation packages if needed
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  BiocManager::install("AnnotationDbi")
}

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}

library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(tibble)

# Map Entrez IDs to gene symbols
feature_annotation <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(entrez_ids),
  keytype = "ENTREZID",
  columns = c("SYMBOL", "GENENAME")
) %>%
  distinct(ENTREZID, .keep_all = TRUE)

head(feature_annotation)

expression_df <- as.data.frame(
  expression_matrix,
  check.names = FALSE
) %>%
  rownames_to_column("Platform_ID") %>%
  mutate(
    ENTREZID = sub("_at$", "", Platform_ID)
  ) %>%
  left_join(
    feature_annotation,
    by = "ENTREZID"
  )

# Check mapping success
table(is.na(expression_df$SYMBOL))
head(expression_df[, c("Platform_ID", "ENTREZID", "SYMBOL", "GENENAME")])

sample_columns <- colnames(expression_matrix)

expression_gene <- expression_df %>%
  dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
  dplyr::select(Gene = SYMBOL, dplyr::all_of(sample_columns)) %>%
  dplyr::mutate(
    dplyr::across(dplyr::all_of(sample_columns), as.numeric)
  ) %>%
  dplyr::group_by(Gene) 


# -------------------------------------------------------------------------
# Extract and clean sample metadata
# -------------------------------------------------------------------------

metadata_raw <- Biobase::pData(gset) %>%
  as.data.frame(check.names = FALSE)

if (!"geo_accession" %in% colnames(metadata_raw)) {
  metadata_raw$geo_accession <- rownames(metadata_raw)
}

# Ensure metadata are in the same order as expression columns.
metadata_raw <- metadata_raw[
  match(colnames(expression_matrix), metadata_raw$geo_accession),
  ,
  drop = FALSE
]

# Identify GEO characteristics columns.
characteristic_columns <- grep(
  "^characteristics_ch1",
  colnames(metadata_raw),
  value = TRUE
)

# -------------------------------------------------------------------------
# Extract and combine metadata from all four GSE36133 ExpressionSet objects
# -------------------------------------------------------------------------

metadata_list <- lapply(gset_list, function(eset) {
  
  metadata_block <- Biobase::pData(eset)
  
  data.frame(
    Sample_ID = as.character(metadata_block$geo_accession),
    Cell_Line = as.character(metadata_block$title),
    Source = as.character(metadata_block$source_name_ch1),
    
    Primary_Site = sub(
      "^primary site:\\s*",
      "",
      as.character(metadata_block$characteristics_ch1)
    ),
    
    Histology = sub(
      "^histology:\\s*",
      "",
      as.character(metadata_block$characteristics_ch1.1)
    ),
    
    Histology_Subtype = sub(
      "^histology subtype1:\\s*",
      "",
      as.character(metadata_block$characteristics_ch1.2)
    ),
    
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})

metadata <- do.call(rbind, metadata_list)

# Convert missing or blank subtype values to NA
metadata$Histology_Subtype[
  metadata$Histology_Subtype == ""
] <- NA

# Remove duplicated samples, if any
metadata <- metadata[
  !duplicated(metadata$Sample_ID),
  ,
  drop = FALSE
]

rownames(metadata) <- NULL

# Inspect the result
dim(metadata)
head(metadata)
table(metadata$Primary_Site)

expression_gene
metadata

# ---- perform PHeatmap ----

# ============================================================
# pheatmap clustering using the top 1,000 variable genes
# Inputs:
#   expression_gene: first column = Gene; remaining columns = samples
#   metadata: first column or Sample_ID column = sample identifiers
# ============================================================

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

library(pheatmap)

# ------------------------------------------------------------
# 1. Convert expression_gene to a numeric matrix
# ------------------------------------------------------------

expression_df <- as.data.frame(
  expression_gene,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Assume the first column contains gene symbols
gene_column <- colnames(expression_df)[1]

gene_names <- as.character(expression_df[[gene_column]])

expression_matrix <- as.matrix(
  expression_df[, -1, drop = FALSE]
)

storage.mode(expression_matrix) <- "numeric"
rownames(expression_matrix) <- gene_names

# Remove genes with missing or blank names
keep_gene_name <- !is.na(rownames(expression_matrix)) &
  rownames(expression_matrix) != ""

expression_matrix <- expression_matrix[
  keep_gene_name,
  ,
  drop = FALSE
]

# Remove duplicate genes, retaining the first occurrence.
# expression_gene should already contain one row per gene.
expression_matrix <- expression_matrix[
  !duplicated(rownames(expression_matrix)),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 2. Prepare metadata
# ------------------------------------------------------------

metadata_df <- as.data.frame(
  metadata,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Use Sample_ID when present; otherwise use the first column
if ("Sample_ID" %in% colnames(metadata_df)) {
  sample_id_column <- "Sample_ID"
} else if ("geo_accession" %in% colnames(metadata_df)) {
  sample_id_column <- "geo_accession"
} else {
  sample_id_column <- colnames(metadata_df)[1]
}

metadata_df[[sample_id_column]] <- as.character(
  metadata_df[[sample_id_column]]
)

# Remove duplicate metadata rows
metadata_df <- metadata_df[
  !duplicated(metadata_df[[sample_id_column]]),
  ,
  drop = FALSE
]

rownames(metadata_df) <- metadata_df[[sample_id_column]]

# ------------------------------------------------------------
# 3. Match expression samples and metadata
# ------------------------------------------------------------

common_samples <- intersect(
  colnames(expression_matrix),
  rownames(metadata_df)
)

if (length(common_samples) < 2) {
  stop(
    "Fewer than two matching sample IDs were found between ",
    "expression_gene and metadata."
  )
}

message("Expression samples: ", ncol(expression_matrix))
message("Metadata samples: ", nrow(metadata_df))
message("Matching samples: ", length(common_samples))

# Retain samples in expression-matrix order
expression_matrix <- expression_matrix[
  ,
  common_samples,
  drop = FALSE
]

metadata_df <- metadata_df[
  common_samples,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(expression_matrix),
    rownames(metadata_df)
  )
)

# ------------------------------------------------------------
# 4. Remove unusable genes
# ------------------------------------------------------------

# Retain genes with at least two finite values
keep_finite <- apply(
  expression_matrix,
  1,
  function(x) sum(is.finite(x)) >= 2
)

expression_matrix <- expression_matrix[
  keep_finite,
  ,
  drop = FALSE
]

# Replace non-finite values with NA before calculating variance
expression_matrix[!is.finite(expression_matrix)] <- NA_real_

# Calculate gene variance across samples
gene_variance <- apply(
  expression_matrix,
  1,
  var,
  na.rm = TRUE
)

# Remove genes with zero, missing, or infinite variance
keep_variable <- is.finite(gene_variance) &
  gene_variance > 0

expression_matrix <- expression_matrix[
  keep_variable,
  ,
  drop = FALSE
]

gene_variance <- gene_variance[keep_variable]

# ------------------------------------------------------------
# 5. Select the top 1,000 variable genes
# ------------------------------------------------------------

number_top_genes <- min(
  1000,
  nrow(expression_matrix)
)

top_gene_names <- names(
  sort(
    gene_variance,
    decreasing = TRUE
  )
)[seq_len(number_top_genes)]

heatmap_matrix <- expression_matrix[
  top_gene_names,
  ,
  drop = FALSE
]

message(
  "Using ",
  nrow(heatmap_matrix),
  " most variable genes."
)

# ------------------------------------------------------------
# 6. Handle missing values
#
# Replace any missing expression value with the median expression
# of that gene before scaling and clustering.
# ------------------------------------------------------------

for (i in seq_len(nrow(heatmap_matrix))) {
  
  missing_values <- is.na(heatmap_matrix[i, ])
  
  if (any(missing_values)) {
    heatmap_matrix[i, missing_values] <- median(
      heatmap_matrix[i, ],
      na.rm = TRUE
    )
  }
}

# ------------------------------------------------------------
# 7. Create column annotations
# ------------------------------------------------------------

# Recommended GSE36133 annotation variables
preferred_annotation_columns <- c(
  #"Primary_Site",
  "Histology"
  #"Histology_Subtype"
)

annotation_columns <- intersect(
  preferred_annotation_columns,
  colnames(metadata_df)
)

# If those columns do not exist, use informative categorical columns
if (length(annotation_columns) == 0) {
  
  candidate_columns <- setdiff(
    colnames(metadata_df),
    sample_id_column
  )
  
  informative_columns <- candidate_columns[
    vapply(
      metadata_df[, candidate_columns, drop = FALSE],
      function(x) {
        number_unique <- length(unique(x[!is.na(x) & x != ""]))
        number_unique >= 2 && number_unique <= 50
      },
      logical(1)
    )
  ]
  
  annotation_columns <- head(informative_columns, 3)
}

if (length(annotation_columns) > 0) {
  
  annotation_col <- metadata_df[
    ,
    annotation_columns,
    drop = FALSE
  ]
  
  # pheatmap handles categorical annotations most reliably as factors
  annotation_col[] <- lapply(
    annotation_col,
    function(x) {
      x <- as.character(x)
      x[is.na(x) | x == ""] <- "Unknown"
      factor(x)
    }
  )
  
} else {
  
  annotation_col <- NULL
  warning("No suitable metadata annotation columns were identified.")
}

# ------------------------------------------------------------
# 8. Row Z-score normalization
#
# pheatmap(scale = "row") centers and scales each gene across
# all samples, allowing clustering by relative expression pattern.
# ------------------------------------------------------------
pheatmap::pheatmap(
  heatmap_matrix,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  show_rownames = FALSE,
  show_colnames = FALSE
)



# ============================================================
# ----UMAP of GSE36133 using the top 1000 variable genes ----
# ============================================================

if (!requireNamespace("uwot", quietly = TRUE))
  install.packages("uwot")

if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

library(uwot)
library(ggplot2)

# ------------------------------------------------------------
# Expression matrix
# ------------------------------------------------------------

expr <- expression_gene

gene_names

# Extract the first column as a vector
gene_names <- as.character(expression_df[[1]])

# Convert all remaining columns to a numeric matrix
expr_mat <- as.matrix(
  expression_df[, -1, drop = FALSE]
)

storage.mode(expr_mat) <- "numeric"

rownames(expr_mat) <- gene_names

# ------------------------------------------------------------
# Remove duplicated genes
# ------------------------------------------------------------

expr_mat <- expr_mat[
  !duplicated(rownames(expr_mat)),
]

# ------------------------------------------------------------
# Select top 1000 variable genes
# ------------------------------------------------------------

gene_var <- apply(expr_mat,1,var,na.rm=TRUE)

top1000 <- names(sort(gene_var,decreasing=TRUE))[1:min(1000,length(gene_var))]

expr_mat <- expr_mat[top1000,]

# ------------------------------------------------------------
# Scale each gene
# ------------------------------------------------------------

expr_scaled <- t(scale(t(expr_mat)))

expr_scaled[is.na(expr_scaled)] <- 0

# ------------------------------------------------------------
# PCA
# ------------------------------------------------------------

pca <- prcomp(
  t(expr_scaled),
  center=FALSE,
  scale.=FALSE
)

# ------------------------------------------------------------
# UMAP
# ------------------------------------------------------------

set.seed(123)

umap <- uwot::umap(
  pca$x[,1:20],
  n_neighbors=30,
  min_dist=0.3,
  metric="euclidean"
)

umap <- as.data.frame(umap)

colnames(umap) <- c("UMAP1","UMAP2")

# ------------------------------------------------------------
# Merge metadata
# ------------------------------------------------------------

sample_col <- intersect(
  c("Sample_ID","geo_accession"),
  colnames(metadata)
)[1]

umap$Sample_ID <- rownames(pca$x)

plot_df <- merge(
  umap,
  metadata,
  by.x="Sample_ID",
  by.y=sample_col,
  all.x=TRUE,
  sort=FALSE
)

# ------------------------------------------------------------
# Plot by histology
# ------------------------------------------------------------



library(ggplot2)
library(ggrepel)
library(dplyr)

# Keep only valid UMAP coordinates
plot_df_clean <- plot_df %>%
  dplyr::mutate(
    UMAP1 = as.numeric(UMAP1),
    UMAP2 = as.numeric(UMAP2),
    Histology = as.character(Histology)
  ) %>%
  dplyr::filter(
    is.finite(UMAP1),
    is.finite(UMAP2),
    !is.na(Histology),
    Histology != ""
  )

# Calculate one label position per histology
histology_labels <- plot_df_clean %>%
  dplyr::group_by(Histology) %>%
  dplyr::summarise(
    UMAP1 = median(UMAP1, na.rm = TRUE),
    UMAP2 = median(UMAP2, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(
    is.finite(UMAP1),
    is.finite(UMAP2),
    n >= 3
  )



umap_basic <- ggplot(
  plot_df_clean,
  aes(
    x = UMAP1,
    y = UMAP2,
    color = Histology
  )
) +
  geom_point(size = 2, alpha = 0.8) +
  theme_bw(base_size = 14) +
  labs(
    title = "GSE36133 UMAP",
    subtitle = "Top 1,000 variable genes",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Histology"
  )

umap_labeled <- umap_basic +
  geom_text(
    data = histology_labels,
    mapping = aes(
      x = UMAP1,
      y = UMAP2,
      label = Histology
    ),
    inherit.aes = FALSE,
    color = "black",
    size = 4,
    fontface = "bold",
    check_overlap = TRUE
  )
umap_labeled
