rm(list = ls())
library("edgeR")
library("limma")
library("Glimma")
library("gplots")
library("org.Hs.eg.db")
library("RColorBrewer")
library("ggplot2")

########## Create heatmap of ARTarget-BMIL genes
TEDNrawdata2 <- read.delim("team4/data/ARTarget-BMIL.csv",
    check.names = FALSE, stringsAsFactors = FALSE
)
dim(TEDNrawdata2)

TEDN2 <- DGEList(
    counts = TEDNrawdata2[, 2:5],
    genes = TEDNrawdata2[, 1]
)
dim(TEDN2)

logcounts <- cpm(TEDN2, log = TRUE)
var_genes <- apply(logcounts, 1, var)

select_var <- names(var_genes)

highly_variable_lcpm <- logcounts[select_var, ]
dim(highly_variable_lcpm)

mypalette <- brewer.pal(11, "PuOr")
morecols <- colorRampPalette(mypalette)

pdf(file = "team4/results/ARTarget-BMIL.heatmap.pdf")
heatmap.2(highly_variable_lcpm,
    trace = "none",
    col = rev(morecols(400)), main = "",
    scale = "row",
    Colv = "NA", Rowv = "NA", dendrogram = "none", cexCol = 1.2,
    density.info = "none", key.title = "NA", key.xlab = " ",
    lmat = rbind(c(0, 3, 4), c(2, 1, 0)), lwid = c(0.5, 4, 1),
    lhei = c(1, 4), keysize = 0.5
)
dev.off()


########## Create heatmap of ARTarget-RNF2 genes
TEDNrawdata3 <- read.delim("team4/data/ARTarget-RNF2.csv",
    check.names = FALSE, stringsAsFactors = FALSE
)
dim(TEDNrawdata3)

TEDN3 <- DGEList(
    counts = TEDNrawdata3[, 2:5],
    genes = TEDNrawdata3[, 1]
)
dim(TEDN3)

logcounts <- cpm(TEDN3, log = TRUE)
var_genes <- apply(logcounts, 1, var)
head(var_genes)

select_var <- names(var_genes)
head(select_var)
highly_variable_lcpm <- logcounts[select_var, ]
dim(highly_variable_lcpm)

mypalette <- brewer.pal(11, "PuOr")
morecols <- colorRampPalette(mypalette)

pdf(file = "team4/results/ARTarget-RNF2.heatmap.pdf")
heatmap.2(highly_variable_lcpm,
    trace = "none",
    col = rev(morecols(400)), main = "",
    scale = "row",
    Colv = "NA", Rowv = "NA", dendrogram = "none", cexCol = 1.2,
    density.info = "none", key.title = "NA", key.xlab = " ",
    lmat = rbind(c(0, 3, 4), c(2, 1, 0)), lwid = c(0.5, 4, 1),
    lhei = c(1, 4), keysize = 0.5
)
dev.off()
