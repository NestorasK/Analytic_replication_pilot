##############################################################################
# The data comes from a Nature Communications paper,
# "BMI1 regulates androgen receptor in prostate cancer independently of the
# polycomb repressive complex 1" (Zhu et al. 2018). Both the raw data
# (sequence reads) and processed data (counts) can be downloaded from GEO under
# accession number GSE97831.

# This study examines the expression profiles of BMI1 knockdown or RING1B
# knockdown in prostate cancer cell line C4-2.Three groups are present,
# and each group contains two biological replicates. I used the counts file
# (tab-delimited text file contian normalized count values for each gene of
# each sample. CPM (counts per million) values) above 0.5) as a starting point
# for analysis. This data has already been aligned to the human gnome.
# https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html




########## Load all packages I need to analysis the data.
library("edgeR")
library("limma")
library("Glimma")
library("gplots")
library("org.Hs.eg.db")
library("RColorBrewer")
library("ggplot2")

########## --------Reading in the data
rawdata <- read.delim("team4/data/GSE97831_genes_expression_table.txt",
    check.names = FALSE, stringsAsFactors = FALSE
)

######## -------For easy manipulation, convert counts data to a DGEList object:
y <- DGEList(counts = rawdata[, 2:7], genes = rawdata[, 1])

########## -------Quality control
y$samples$lib.size
barplot(y$samples$lib.size,
    main = "Barplot of library sizes",
    names = colnames(y), las = 2
)

########## -------Distribution plots
logcounts <- cpm(y, log = TRUE)
boxplot(logcounts,
    xlab = "", ylab = "Log2 counts per million",
    las = 2
)
abline(h = median(logcounts), col = "blue")
title("Boxplots of logCPMs (unnormalised)")

########## --------Multidimensional scaling plots
col.cell <- c(
    "blue", "blue", "yellow", "yellow",
    "Green", "Green"
)
plotMDS(y, col = col.cell)

########## -------Hierarchical clustering with heatmaps
var_genes <- apply(logcounts, 1, var)
select_var <- names(sort(var_genes, decreasing = TRUE))[1:500]
head(select_var)
highly_variable_lcpm <- logcounts[select_var, ]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

mypalette <- brewer.pal(11, "RdYlBu")
morecols <- colorRampPalette(mypalette)

png(file = "team4/results/High_var_genes500.heatmap.png")
heatmap.2(highly_variable_lcpm,
    trace = "none",
    col = rev(morecols(50)),
    main = "Top 500 most variable genes",
    ColSideColors = col.cell, scale = "row"
)
dev.off()

########## --------Normalisation for composition bias
y <- calcNormFactors(y)
y$samples
par(mfrow = c(1, 2))
plotMD(logcounts, column = 3)
abline(h = 0, col = "grey")
plotMD(logcounts, column = 4)
abline(h = 0, col = "grey")
plotMD(logcounts, column = 6)
abline(h = 0, col = "grey")

########## ---------Differential expression
group <- factor(c(
    "siCont", "siCont", "siBMIL",
    "siBMIL", "siRNF2", "siRNF2"
), levels = c("siCont", "siBMIL", "siRNF2"))
design <- model.matrix(~ 0 + group)

par(mfrow = c(1, 1))
v <- voom(y, design, plot = TRUE)

fit <- lmFit(v)

cont.matrix <- makeContrasts(
    BMILVsCont = groupsiBMIL - groupsiCont,
    RNF2VsCont = groupsiRNF2 - groupsiCont,
    levels = design
)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
summa.fit <- decideTests(fit.cont)
summary(summa.fit)
results <- decideTests(fit.cont)
a <- vennCounts(results)
print(a)
mfrow.old <- par()$mfrow
par(mfrow = c(1, 2))
vennDiagram(a)
vennDiagram(results,
    include = c("up", "down"),
    counts.col = c("red", "blue"),
    circle.col = c("red", "blue", "green3")
)
par(mfrow = mfrow.old)
topTable(fit.cont, coef = "BMILVsCont", sort.by = "p")
topTable(fit.cont, coef = "RNF2VsCont", sort.by = "p")

########## Adding annotation and saving the results
limma.res <- topTable(fit.cont, coef = "BMILVsCont", sort.by = "p", n = "Inf")
write.csv(limma.res, file = "team4/results/BMILVsCon.csv", row.names = FALSE)
limma.res <- topTable(fit.cont, coef = "RNF2VsCont", sort.by = "p", n = "Inf")
write.csv(limma.res, file = "team4/results/RNF2VsCont.csv", row.names = FALSE)

########## ------Plots after testing for DE
volcanoplot(fit.cont, coef = 1, highlight = 100, names = fit.cont$genes$genes)
volcanoplot(fit.cont, coef = 2, highlight = 100, names = fit.cont$genes$genes)

########## ------Create heatmap of DE genes
TEDNrawdata <- read.delim("team4/data/TE.txt",
    check.names = FALSE,
    stringsAsFactors = FALSE
)
dim(TEDNrawdata)
TEDN <- DGEList(counts = TEDNrawdata[, 2:7], genes = TEDNrawdata[, 1])
dim(TEDN)

logcounts <- cpm(TEDN, log = TRUE)
var_genes <- apply(logcounts, 1, var)
select_var <- names(var_genes)
highly_variable_lcpm <- logcounts[select_var, ]
dim(highly_variable_lcpm)

mypalette <- brewer.pal(11, "PuOr")
morecols <- colorRampPalette(mypalette)

png(file = "team4/results/TEDN-F.heatmap.png")
heatmap.2(highly_variable_lcpm,
    trace = "none",
    col = rev(morecols(400)),
    main = "",
    scale = "row",
    Colv = "NA",
    Rowv = "NA",
    dendrogram = "none",
    cexCol = 1.2,
    density.info = "none",
    key.title = "NA",
    key.xlab = " ",
    lmat = rbind(c(0, 3, 4), c(2, 1, 0)),
    lwid = c(0.5, 4, 1),
    lhei = c(1, 4),
    keysize = 0.5
)
dev.off()
