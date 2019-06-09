# CHETAH paper figure scripts
This repository contains all the scripts necessary for producing the figures of the CHETAH article, for which a preprint is available at: https://www.biorxiv.org/content/10.1101/558908v1

## versions
CHETAH is currently available at Bioconductor: https://www.bioconductor.org/packages/release/bioc/html/CHETAH.html
These scripts were produced with an earlier version of CHETAH, which uses exactly the same code as the Bioconductor version, but did not yet integrate the compatibility with the “SingleCellExperiment” package. They were run in R version 3.5.1 and Seurat version 2.3.4 was used 
This version of CHETAH and Seurat can be installed by running: 

```{r echo=TRUE, eval=FALSE}
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("jdekanter/CHETAH",  ref = “c0091a1”)

## Check if all dependencies are installed
dep <- c('bioDist', 'ggplot2', 'gplots', 'cowplot',
         'dendextend', 'corrplot', 'reshape2', 'plotly', 'grDevices')
pkg_avail <- suppressMessages(sapply(dep, function (pkg) pkg %in% installed.packages()[, "Package"]))

# --- Install dependencies, if neccesary
if(length(dep[!pkg_avail]) > 0) {
  if (!require("BiocManager")) {
    install.packages("BiocManager")
  }
  BiocManager::install(dep[!pkg_avail])
}
# Load the package
library(CHETAH)
## Install Seurat
devtools::install_version(package = 'Seurat', version = package_version('2.3.4'))
```

## Running the scripts
First run all the “Import_*.R” scripts. The location of the raw data is mentioned in each script. Also download the Tumor_reference from: https://figshare.com/s/aaf026376912366f81b6 
Then run all the “Analysis_compare.R” script: this is necessary for Figure 2, 3C and S2-5.
Finally, run all other “Figure*.R” scripts.
