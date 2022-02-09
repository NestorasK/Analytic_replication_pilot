# Figure 2E - NOT REPRODUCED #####
# Data available in file: 1-s2.0-S0092867411007677-mmc2.xls, sheet B
# of the initial publication
# Export sheet B to the Figure2E.csv file
rm(list = ls())
library("gplots")
setwd("PATH_TO_Analytic_replication_pilot")

spread <- read.csv(
    file = "team1/data/Figure2E.csv", header = TRUE,
    stringsAsFactors = TRUE
)
row.names(spread) <- spread$gene
spread <- spread[, 2:17]
datamatrix <- data.matrix(spread)
tdatamatrix <- t(datamatrix)

color <- colorRampPalette(c("black", "red"))
png(
    filename = "team1/results/Figure2E_1.png", width = 700, height = 700,
    res = 100
)
heatmap.2(datamatrix,
    Rowv = TRUE, Colv = TRUE, scale = "column",
    dendrogram = "both", col = color, key = TRUE,
    density.info = "none", trace = "none", key.title = "",
    keysize = 1, key.xlab = "log rpkmx10^3", margins = c(12, 8),
    cexCol = 1.2
)
dev.off()

# Determining Pearson's correlation for clustering pattern:
# Determining Pearson's correlation between tissues (columns)
r2tissue <- cor(
    x = datamatrix, use = "pairwise.complete.obs",
    method = "pearson"
)
tissueclust <- hclust(dist(r2tissue), method = "complete")
plot(tissueclust, main = "", xlab = "Tissue")

# Determining Pearson's correlation between genes (rows)
r2gene <- cor(
    x = tdatamatrix, use = "pairwise.complete.obs",
    method = "pearson"
)
geneclust <- hclust(dist(r2gene), method = "complete")
plot(geneclust, main = "", xlab = "Gene")

# Making a heatmap using the Pearson's correlation for clustering
png(
    filename = "team1/results/Figure2E_using_Pearson.png",
    width = 700, height = 700,
    res = 100
)
heatmap.2(datamatrix,
    Rowv = TRUE, Colv = TRUE, scale = "column",
    dendrogram = "both", col = color, key = TRUE, density.info = "none",
    trace = "none", key.title = "", keysize = 1,
    key.xlab = "log rpkmx10^3",
    distfun = function(datamatrix) as.dist(1 - cor(t(datamatrix))),
    hclustfun = function(datamatrix) hclust(datamatrix, method = "complete"),
    margins = c(12, 8),
    cexCol = 1.2
)
dev.off()
