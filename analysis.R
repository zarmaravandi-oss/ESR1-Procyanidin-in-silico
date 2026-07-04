# ============================================================
# Transcriptomic analysis workflow related to:
# "In silico evaluation of procyanidin as a potential ESR1 inhibitor:
# docking and MD insights in uterine fibroids and endometriosis"
#
# This script contains the retained R analysis steps used in the project,
# cleaned for readability and repository sharing.
#
# Notes:
# - Raw GEO data were downloaded previously.
# - Some normalization / preprocessing steps were performed earlier and are
#   not fully preserved in this script.
# - This file reflects the practical downstream workflow retained from the
#   transcriptomic analyses used in the study.
#
# Main components:
# - ESR1 expression analysis
# - statistical comparison across groups
# - GEO retrieval for fibroid dataset (GSE64763)
# - probe-to-gene mapping
# - DEG-based enrichment analysis (GO / KEGG)
# - export of supplementary enrichment tables
# ============================================================


# ============================================================
# 1) Package installation (if needed)
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c(
  "dplyr",
  "readr",
  "ggplot2"
)

bioc_packages <- c(
  "GEOquery",
  "limma",
  "affy",
  "oligo",
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "hgu133plus2.db"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}


# ============================================================
# 2) Load libraries
# ============================================================

library(GEOquery)
library(limma)
library(affy)
library(oligo)

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(hgu133plus2.db)

library(dplyr)
library(readr)
library(ggplot2)


# ============================================================
# 3) ESR1 expression analysis in endometriosis dataset
# ============================================================

esr1_expr <- expr_mat["ESR1", ]

group <- meta_table$Group
names(group) <- meta_table$SampleName

tapply(esr1_expr, group, mean)

t.test(
  esr1_expr[group == "DiEIn"],
  esr1_expr[group == "PE"]
)


# ============================================================
# 4) ESR1 expression comparison in leiomyoma vs normal myometrium
#    using current expr_mat / meta_table objects
# ============================================================

esr1_expr <- expr_mat["ESR1", ]

group <- meta_table$Group
names(group) <- meta_table$SampleName

tapply(esr1_expr, group, mean)

t.test(
  esr1_expr[group == "leiomyoma"],
  esr1_expr[group == "normal myometrium"]
)


# ============================================================
# 5) Download and prepare fibroid GEO dataset: GSE64763
# ============================================================

gse_fib <- getGEO("GSE64763", GSEMatrix = TRUE)
eset_fib <- gse_fib[[1]]

expr_mat_fib <- exprs(eset_fib)
meta_fib <- pData(eset_fib)
fdata_fib <- fData(eset_fib)


# ============================================================
# 6) ESR1 probe mapping and expression analysis in fibroid dataset
# ============================================================

group_fib <- meta_fib$source_name_ch1
names(group_fib) <- rownames(meta_fib)

esr1_idx <- which(fdata_fib$`Gene Symbol` == "ESR1")

esr1_expr_fib_all <- expr_mat_fib[esr1_idx, , drop = FALSE]
esr1_expr_fib <- colMeans(esr1_expr_fib_all)

fib_vs_norm <- t.test(
  esr1_expr_fib[group_fib == "Fibroid"],
  esr1_expr_fib[group_fib == "Normal myometrium"],
  var.equal = FALSE
)

fib_vs_norm


# ============================================================
# 7) DEG object preparation
# ============================================================

degs_endo <- degs_PE
degs_endo$Gene.symbol <- rownames(degs_PE)

deg_ulms_vs_normal$ID <- rownames(deg_ulms_vs_normal)

probe_ids <- rownames(deg_ulms_vs_normal)

gene_symbols <- mapIds(
  hgu133plus2.db,
  keys = probe_ids,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

deg_ulms_vs_normal$GeneSymbol <- gene_symbols

degs_endo_sig <- degs_endo %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

deg_fib <- deg_ulms_vs_normal
deg_fib_sig <- deg_fib %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 0.5)


# ============================================================
# 8) Entrez ID conversion
# ============================================================

genes_endo <- bitr(
  degs_endo_sig$Gene.symbol,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

genes_fib <- bitr(
  deg_fib_sig$GeneSymbol,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)


# ============================================================
# 9) KEGG enrichment
# ============================================================

kegg_endo <- enrichKEGG(
  gene = genes_endo$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

kegg_fib <- enrichKEGG(
  gene = genes_fib$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)


# ============================================================
# 10) GO enrichment
# ============================================================

go_endo <- enrichGO(
  gene = genes_endo$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

go_fib <- enrichGO(
  gene = genes_fib$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE,
  qvalueCutoff = 0.25
)


# ============================================================
# 11) Enrichment plots
# ============================================================

barplot(go_endo, showCategory = 20, title = "GO BP - Endometriosis")
dotplot(go_fib, showCategory = 20, title = "GO BP - Fibroids")


# ============================================================
# 12) Export supplementary enrichment tables
# ============================================================

outdir <- "supplementary_enrichment"
if (!dir.exists(outdir)) dir.create(outdir)

if (nrow(as.data.frame(kegg_endo)) > 0) {
  write_csv(
    as.data.frame(kegg_endo),
    file.path(outdir, "SuppTable_KEGG_Endometriosis.csv")
  )
}

if (nrow(as.data.frame(go_endo)) > 0) {
  write_csv(
    as.data.frame(go_endo),
    file.path(outdir, "SuppTable_GO_Endometriosis.csv")
  )
}

if (nrow(as.data.frame(kegg_fib)) > 0) {
  write_csv(
    as.data.frame(kegg_fib),
    file.path(outdir, "SuppTable_KEGG_Fibroid.csv")
  )
}

if (nrow(as.data.frame(go_fib)) > 0) {
  write_csv(
    as.data.frame(go_fib),
    file.path(outdir, "SuppTable_GO_Fibroid.csv")
  )
}

message("Analysis completed. CSV files were saved in 'supplementary_enrichment'.")

