# Installing R and RStudio

Before running this pipeline, install **R** and **RStudio**.

## Step 1. Install R

Download the latest version of R from the Comprehensive R Archive Network (CRAN):

https://cran.r-project.org/

Choose your operating system:

* **Windows:** Download and run the Windows installer.
* **macOS:** Download the macOS installer (.pkg) and complete the installation.
* **Linux:** Follow the installation instructions for your distribution on the CRAN website.

After installation, verify that R is installed by opening a terminal (or Command Prompt) and running:

```bash
R --version
```

---

## Step 2. Install RStudio

Download the free **RStudio Desktop** from Posit:

https://posit.co/download/rstudio-desktop/

Select the installer for your operating system and follow the installation instructions.

After installation:

1. Open **RStudio**.
2. The R Console should appear automatically.
3. Confirm the R version by running:

```r
R.version.string
```

---

## Step 3. Install Required Packages

Open RStudio and execute:

```r
install.packages(c(
  "dplyr",
  "tidyr",
  "tibble",
  "readr",
  "stringr",
  "ggplot2",
  "ggrepel",
  "uwot",
  "pheatmap",
  "BiocManager"
))
```

Next, install the required Bioconductor packages:

```r
BiocManager::install(c(
  "GEOquery",
  "Biobase",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "biomaRt"
))
```

---

## Step 4. Verify the Installation

Run the following commands in the R Console:

```r
library(GEOquery)
library(dplyr)
library(ggplot2)
library(uwot)
library(pheatmap)
```

If no errors are reported, your installation is complete.

---

## Step 5. Running the Script

Open the R script in RStudio and click **Source**, or run it from a terminal:

```bash
Rscript download_GSE36133.R
```

The output files will be written to the `GSE36133_Output/` directory.

---

## Troubleshooting

### Package installation fails

Update R to the latest version and try again.

### Bioconductor package not found

Update Bioconductor:

```r
BiocManager::install(version = "3.22")
BiocManager::install()
```

### Unable to download GEO datasets

Ensure your internet connection is active and that the GEO servers are accessible.

### Need help?

Visit:

* CRAN: https://cran.r-project.org/
* Posit Community: https://forum.posit.co/
* Bioconductor: https://bioconductor.org/
