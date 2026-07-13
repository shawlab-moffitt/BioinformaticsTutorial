# GSE36133 Downloader and Visualization Pipeline

## Overview

This R script downloads the **GSE36133 (Cancer Cell Line Encyclopedia, CCLE)** microarray dataset from GEO and generates analysis-ready expression and metadata tables. It further annotates genes with HGNC symbols, filters to protein-coding genes, performs unsupervised hierarchical clustering, and generates a UMAP representation of the cell lines colored by histology.

The workflow was developed as an example preprocessing pipeline for **Easy-App / DRPPM-EASY**.

---

## Features

* Downloads **GSE36133** directly from GEO using `GEOquery`
* Maps BrainArray Entrez identifiers to HGNC gene symbols
* Produces both:

  * all annotated genes
  * protein-coding genes only
* Extracts standardized sample metadata
* Generates:

  * expression matrices
  * metadata table
  * hierarchical clustering heatmap
  * UMAP visualization

---

## Requirements

### R packages

Bioconductor

* GEOquery
* Biobase
* AnnotationDbi
* org.Hs.eg.db
* biomaRt

CRAN

* dplyr
* tidyr
* tibble
* readr
* stringr
* pheatmap
* uwot
* ggplot2
* ggrepel

---

## Output Files

The script creates the directory

```
GSE36133_Output/
```

and writes

| File                                            | Description                                                |
| ----------------------------------------------- | ---------------------------------------------------------- |
| `GSE36133_gene_expression.tsv`                  | Gene-level expression matrix mapped to HGNC symbols        |
| `GSE36133_proteincoding_expression.tsv`         | Protein-coding gene expression matrix                      |
| `GSE36133_meta.tsv`                             | Sample metadata                                            |
| `GSE36133_top1000_variable_genes_pheatmap.pdf`* | Hierarchical clustering heatmap (if PDF output is enabled) |
| `GSE36133_UMAP_Histology_labeled.png/pdf`*      | UMAP visualization (optional)                              |

* Generated when the corresponding save commands are enabled.

---

## Metadata Fields

The metadata table contains

* Sample_ID
* Cell_Line
* Source
* Primary_Site
* Histology
* Histology_Subtype

---

## Expression Processing

The pipeline

1. Downloads the GEO Series Matrix.
2. Maps Entrez identifiers to HGNC gene symbols.
3. Removes genes without valid symbols.
4. Collapses duplicated symbols by averaging expression values.
5. Retrieves protein-coding genes from Ensembl using **biomaRt**.
6. Produces a filtered protein-coding expression matrix.

---

## Heatmap Workflow

The heatmap is generated using:

* Top 1,000 most variable genes
* Row Z-score normalization
* Euclidean distance
* Complete-linkage hierarchical clustering
* Sample annotations from metadata

Column annotations currently use:

* Histology

but can easily be changed to

* Primary_Site
* Histology_Subtype
* or any metadata column.

---

## UMAP Workflow

The dimensionality reduction pipeline consists of

1. Top 1,000 variable genes
2. Gene-wise Z-score scaling
3. Principal component analysis
4. UMAP computed from the first 20 principal components

The resulting UMAP is colored by histology and can optionally display histology labels.

---

## Running

```bash
Rscript download_GSE36133.R
```

---

## Dataset

**GSE36133**

Barretina et al.

Cancer Cell Line Encyclopedia (CCLE)

Approximately 900 cancer cell lines profiled on the BrainArray GPL15308 platform.

---

## Citation

If this script contributes to your work, please cite:

* GEO: GSE36133
* Cancer Cell Line Encyclopedia (CCLE)
* Easy-App / DRPPM-EASY
