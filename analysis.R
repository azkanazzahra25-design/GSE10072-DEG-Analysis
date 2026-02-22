x <- 10 
x + 5
x > 5
#1. Cek apakah BiocManager sudah ada, jika belum instal dari CRAN
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#2. Instal versi dasar Bioconductor
BiocManager::install(version = "3.22") # Sesuaikan versi dengan R Anda
library(BiocManager)
BiocManager::install(c("GEOquery","limma"))
install.packages(c("pheatmap","ggplot2","dplyr"))
library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(gplots)
library

install.packages("gplots")
library(illuminaHumanv4.db)
library(AnnotationDbi)
library(umap)

BiocManager::install(c("illuminaHumanv4.db","AnnotationDbi"))
install.packages("umap")
gset <- getGEO("GSE10072", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm =
                            TRUE))
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm =
                            TRUE))
options(timeout = 600)
gset <- getGEO("GSE10072", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm =
                            TRUE))
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
ex <- na.omit(ex)
#Jika LogTransform = TRUE, maka lakukan log2
if (LogTransform) {
  # Nilai <= 0 tidak boleh di-log, maka diubah menjadi NA
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}
group_info <- pData(gset)[["source_name_ch1"]]
groups <- make.names(group_info)
gset$group <- factor(groups)

nama_grup <- levels(gset$group)
print(nama_grup)

design <- model.matrix(~0 + gset$group)
colnames(design) <- levels(gset$group)

grup_kanker <- nama_grup[1]
grup_normal <- nama_grup[2]

contrast_formula <- paste(grup_kanker, "-", grup_normal)
print(paste("Kontras yang dianalisis:", contrast_formula))

fit <- lmFit(ex, design)

contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

topTableResults <- topTable(fit2, adjust = "fdr", number = Inf)
topTableResults <- topTable(fit2, adjust = "fdr", number = Inf)
library(limma)
fit <- lmFit(ex, design)

contrast_matrix <- makeContrasts(
  contrasts = contrast_formula,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)

fit2 <- eBayes(fit2)

topTableResults <- topTable(fit2, adjust = "fdr", number = Inf)
ls()
grup_kanker <- nama_grup[1]
grup_normal <- nama_grup[2]

contrast_formula <- paste(grup_kanker, "-", grup_normal)

library(limma)

contrast_matrix <- makeContrasts(
  contrasts = contrast_formula,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

topTableResults <- topTable(fit2, adjust = "fdr", number = Inf)
library(limma)
grup_kanker <- nama_grup[1]
grup_normal <- nama_grup[2]
contrast_formula <- paste(grup_kanker, "-", grup_normal)
contrast_formula
nama_grup
grup_kanker <- nama_grup[1]
grup_normal <- nama_grup[2]

contrast_formula <- paste(grup_kanker, "-", grup_normal)
contrast_formula
levels(gset$group)
contrast_formula <- paste(nama_grup[1], "-", nama_grup[2])
contrast_formula
contrast_formula <- paste(nama_grup[1], "-", nama_grup[2])
contrast_formula
contrast_formula <- "Adenocarcinoma.of.the.Lung-Normal.Lung.Tissue"
contrast_formula
contrast_formula <- "Adenocarcinoma.of.the.Lung-Normal.Lung.Tissue"
print(contrast_formula)
library(limma)

contrast_matrix <- makeContrasts(
  contrasts = contrast_formula,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

topTableResults <- topTable(fit2, number = Inf)
library(ggplot2)

topTableResults$threshold <- as.factor(
  ifelse(topTableResults$adj.P.Val < 0.05 & abs(topTableResults$logFC) >= 1,
         ifelse(topTableResults$logFC > 1, "Up", "Down"),
         "NotSig"))

ggplot(topTableResults, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(size = 1.5) +
  theme_minimal()
library(pheatmap)

top50 <- topTableResults[1:50, ]

heatmap_data <- ex[rownames(top50), ]

pheatmap(heatmap_data,
         scale = "row",
         show_rownames = FALSE,
         show_colnames = TRUE)
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)
sig_genes <- rownames(topTableResults[topTableResults$adj.P.Val < 0.05, ])
ego <- enrichGO(
  gene = sig_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP"
)

barplot(ego, showCategory = 10)
ekegg <- enrichKEGG(
  gene = sig_genes,
  organism = "hsa"
)

barplot(ekegg, showCategory = 10)
ekegg <- enrichKEGG(
  gene = sig_genes,
  organism = "hsa"
)

barplot(ekegg, showCategory = 10)
colnames(topTableResults)
BiocManager::install("hgu133a.db")
library(hgu133a.db)
library(AnnotationDbi)
gene_symbols <- mapIds(
  hgu133a.db,
  keys = rownames(topTableResults),
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

topTableResults$SYMBOL <- gene_symbols
gene_symbols <- mapIds(
  hgu133a.db,
  keys = rownames(topTableResults),
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

topTableResults$SYMBOL <- gene_symbols