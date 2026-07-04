# ESR1-Procyanidin-in-silico
Scripts and workflows for the in silico evaluation of procyanidin as a potential ESR1 inhibitor in uterine fibroids and endometriosis.


## 🧬 Transcriptomic Analysis (`analysis.R`)

The `analysis.R` script provides a complete workflow for processing and analyzing transcriptomic data related to ESR1 expression and its regulation by Procyanidin in uterine fibroids and endometriosis.

### Workflow Steps:
1. **Data Preparation:** GEO data retrieval (e.g., GSE64763 for uterine fibroids).
2. **ESR1 Profiling:** Expression profiling of ESR1 and relevant target genes.
3. **DEG Analysis:** Differential Gene Expression (DEG) analysis using the `limma` package.
4. **Functional Enrichment:** GO (Gene Ontology) and KEGG pathway enrichment analysis to identify significant biological pathways.

### Dependencies
Before running the script, ensure you have the required Bioconductor and CRAN packages installed:
```R
# CRAN packages
install.packages(c("tidyverse", "ggplot2", "ggpubr"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "clusterProfiler", "org.Hs.eg.db", "pathview"))



Usage
To execute the transcriptomic analysis, set your R working directory to the project folder and run the script:

R
source("analysis.R")

Note: This script assumes that the required expression matrix (expr_mat), metadata table (meta_table), and other relevant objects are already loaded in your R environment, or that you have configured the file paths within the script accordingly.

