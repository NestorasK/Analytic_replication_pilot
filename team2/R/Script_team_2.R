# In order to run this please place files
# - SraRunInfo.csv
# - gene_annotation.csv
# - FPKMAll Data.csv
# - FPKMaverage.csv and
# - all files contained in SRRrawcounts.zip into your working directory.

# Set working directory to src_analytic_reproduce folder
# setwd("PATH_TO/src_analytic_reproduce/team2/data/")

# Read SRA file infos
sri <- read.csv("SraRunInfo.csv", stringsAsFactors = FALSE)

### NOTE: Lines 8-30 were used to generate files needed to determine
# differential expression. Files were downloaded with lines 10-20 and
# reads were mapped to mm9 mouse reference genome and counted using tools
# on usegalaxy.org. More details below.

### NOTE: NO NEED TO RUN LINES 13-22 AS THE FINAL OUTPUT HAS BEEN PROVIDED ###
# download SRA files
files <- basename(sri$download_path)
for (i in 1:length(files)) download.file(sri$download_path[i], files[i])

# Assure that all the files has been downloaded successfully,
# they're converted to fastq format
stopifnot(all(file.exists(files)))
for (f in files) {
    cmd <- paste("fastq-dump --split-3", f)
    cat(cmd, "\n") # print the current command
    system(cmd) # invoke command
}

# We used http://usegalaxy.org to run the Tuxedo pipeline.
# https://usegalaxy.org/u/dhwangdxl/h/r-data-import
# We groomed the concatenated fastq files using FASTQ groomer (v 1.0.4),
# aligned reads with Bowtie2 and identify splice junctions with
# TopHat, both using TopHat tool (v2.0.14, Galaxy v 0.9)
# with the following changes to default
# parameters: paired end data, mm9 mouse reference genome.
# For expression level estimation, we used TopHat accepted hits
# with Cufflinks (v 2.2.1, Galaxy v 2.2.1.0), with the following changes to
# default parameters: use reference genome mm9; �Yes�
# count hits compatible with reference RNAs
# only; and �Standard Length Correction� only.
# To generate raw read counts from TopHat accepted hits output,
# we used featureCounts (13), a program bundled into the Rsubread package.
# We used the mm9 mouse reference genome and
# selected the following options: GTF.featureType=�exon�,
# GTF.attrType=�gene_id�, isPairedEnd=TRUE.
# We then use EdgeR to calculate fold change differences between samples

# load EdgeR, Statmod library
library(edgeR)
library(statmod)

# gene annotation
annotation <- read.csv("gene_annotation.csv",
    header = TRUE,
    sep = ",", quote = "\"'",
    stringsAsFactors = FALSE
)

# Generate ENTREZID for genes of interests
relgenes <- c(
    "Fcrls", "Gpr34", "Hexb", "Olfml3", "P2ry12",
    "P2ry13", "Sall1", "Sall3", "Siglech", "Slc2a5", "Sparc",
    "Tmem119", "Aif1", "Csf1r", "Cx3cr1", "Itgam", "Trem2",
    "C1qa", "Cd68", "Il1a", "Tnf"
)
genes <- rep(x = "NA", times = 21)
for (i in 1:length(genes)) {
    genes[i] <- annotation[which(annotation$SYMBOL == relgenes[i]), 1]
}

# read raw counts files generated with Galaxy
files <- dir(pattern = "*\\.txt$")
RG <- readDGE(files) # use readDGE to generate data frame of raw counts
colnames(RG$counts) <- paste(sri$SampleName)
rownames(RG$samples) <- paste(sri$SampleName)

# Generate condition factors
group <- c(
    "ICTBM", "ICTBM", "WTMG", "Cult11B", "BMTKO",
    "Cult11B", "Cult11B", "ICTBlood", "ICTBlood", "ICTBM",
    "ICTBlood", "ICTFetLiv", "ICTFetLiv", "ICTFetLiv", "ICTMG", "ICTMG",
    "ICTFetLiv", "ICTFetLiv", "ICTFetLiv", "ICTYoungMG", "ICTYoungMG",
    "ICTFetBrain", "ICTFetLiv", "ICTFetBrain", "ICTFetBrain",
    "ICTFetBrain", "ICTFetBrain", "ICTCultMG", "ICTCultMG", "ICTCultMG",
    "ICTCultMG", "WTMG", "WTMG", "ICTYS", "WTMG", "ICTYS", "ICTYS",
    "ICTYoungMG", "ICTYoungMG", "WTMG", "WTMG", "BMTKO", "BMTKO", "WTMG",
    "BMTCtl", "Acute11B", "Acute11B", "Acute11B"
)
group1 <- c(group, "BMTKO", "BMTKO", "BMTCtl", "BMTCtl")
RG$samples$group <- as.factor(group1)

# filtering lowly expressed genes with less than 25 counts per million
# in 2 or more samples
cpmReq <- 25
sampleReq <- 2
keep2 <- rowSums(cpm(RG$counts) >= cpmReq) >= sampleReq
RG1 <- RG[keep2, , keep.lib.sizes = FALSE]

# recalculate library size
newlibsize <- rep(NA, times = ncol(RG1$counts))
for (i in 1:ncol(RG1$counts)) {
    newlibsize[i] <- sum(RG1$counts[, i])
}
RG1$samples$lib.size <- newlibsize

# normalization
RG1 <- calcNormFactors(RG1)

# Calculate Dispersion
group1 <- as.factor(group1)
design <- model.matrix(~ 0 + group1, data = RG1$samples)
colnames(design) <- levels(RG1$samples$group)
RG1 <- estimateDisp(RG1, design, robust = TRUE)

# Differential Gene Expression
fit <- glmQLFit(RG1, design, robust = TRUE)
qlf.ICTCultMGvsWTMG <- glmTreat(fit,
    contrast = c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1),
    lfc = log2(2)
)
qlf.ICTMGvsWTMG <- glmTreat(fit,
    contrast = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1),
    lfc = log2(2)
)
qlf.ICTYoungMGvsWTMG <- glmTreat(fit,
    contrast = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1),
    lfc = log2(2)
)
qlf.CultMGvsWTMG <- glmTreat(fit,
    contrast = c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1),
    lfc = log2(1.5)
)

# Testing for differential genes expression for genes of interest for each
# condition
# Generate gene expression data for ICT Cultured Microglia compared to WT
# Microglia
topexpressedgenes <- topTags(qlf.ICTCultMGvsWTMG,
    n = nrow(qlf.ICTCultMGvsWTMG$table)
)
row_names <- rownames(topTags(qlf.ICTCultMGvsWTMG,
    n = nrow(qlf.ICTCultMGvsWTMG$table)
))
dfrICTCultMGvsWTMG <- topexpressedgenes[row_names == genes[1] |
    row_names == genes[2] | row_names == genes[3] | row_names == genes[4] |
    row_names == genes[5] | row_names == genes[6] | row_names == genes[7] |
    row_names == genes[8] | row_names == genes[9] | row_names == genes[10] |
    row_names == genes[11] | row_names == genes[12] | row_names == genes[13] |
    row_names == genes[14] | row_names == genes[15] | row_names == genes[16] |
    row_names == genes[17] | row_names == genes[18] | row_names == genes[19] |
    row_names == genes[20] | row_names == genes[21], ]
ICTCultMGvsWTMGmatrix <- as.matrix(dfrICTCultMGvsWTMG$table)
ICTCultMGvsWTMG <- matrix(rep(NA, times = 21 * 5), nrow = 21, ncol = 5)
for (i in 1:length(genes)) {
    for (j in 1:nrow(ICTCultMGvsWTMGmatrix)) {
        if (genes[i] == rownames(ICTCultMGvsWTMGmatrix)[j]) {
            ICTCultMGvsWTMG[i, ] <- ICTCultMGvsWTMGmatrix[j, ]
        }
    }
}
rownames(ICTCultMGvsWTMG) <- relgenes
colnames(ICTCultMGvsWTMG) <- colnames(ICTCultMGvsWTMGmatrix)

# Generate gene expression data for ICT Adult Microglia compared to WT Microglia
row_names <- rownames(topTags(qlf.ICTMGvsWTMG, n = nrow(qlf.ICTMGvsWTMG$table)))
dfrICTMGvsWTMG <- topTags(qlf.ICTMGvsWTMG,
    n = nrow(qlf.ICTMGvsWTMG$table)
)[row_names == genes[1] | row_names == genes[2] | row_names == genes[3] |
    row_names == genes[4] | row_names == genes[5] | row_names == genes[6] |
    row_names == genes[7] | row_names == genes[8] | row_names == genes[9] |
    row_names == genes[10] | row_names == genes[11] | row_names == genes[12] |
    row_names == genes[13] | row_names == genes[14] | row_names == genes[15] |
    row_names == genes[16] | row_names == genes[17] | row_names == genes[18] |
    row_names == genes[19] | row_names == genes[20] | row_names == genes[21], ]
ICTMGvsWTMGmatrix <- as.matrix(dfrICTMGvsWTMG$table)
ICTMGvsWTMG <- matrix(rep(NA, times = 21 * 5), nrow = 21, ncol = 5)
for (i in 1:length(genes)) {
    for (j in 1:nrow(ICTMGvsWTMGmatrix)) {
        if (genes[i] == rownames(ICTMGvsWTMGmatrix)[j]) {
            ICTMGvsWTMG[i, ] <- ICTMGvsWTMGmatrix[j, ]
        }
    }
}
rownames(ICTMGvsWTMG) <- relgenes
colnames(ICTMGvsWTMG) <- colnames(ICTMGvsWTMGmatrix)

# Generate gene expression data for ICT Young Microglia compared to WT Microglia
topexpressedgenes1 <- topTags(qlf.ICTYoungMGvsWTMG,
    n = nrow(qlf.ICTYoungMGvsWTMG$table)
)
row_names <- rownames(topexpressedgenes1)
dfrICTYoungMGvsWTMG <- topexpressedgenes1[row_names == genes[2] |
    row_names == genes[1] | row_names == genes[3] | row_names == genes[4] |
    row_names == genes[5] | row_names == genes[6] | row_names == genes[7] |
    row_names == genes[8] | row_names == genes[9] | row_names == genes[10] |
    row_names == genes[11] | row_names == genes[12] | row_names == genes[13] |
    row_names == genes[14] | row_names == genes[15] | row_names == genes[16] |
    row_names == genes[17] | row_names == genes[18] | row_names == genes[19] |
    row_names == genes[20] | row_names == genes[21], ]
ICTYoungMGvsWTMGmatrix <- as.matrix(dfrICTYoungMGvsWTMG$table)
ICTYoungMGvsWTMG <- matrix(rep(NA, times = 21 * 5), nrow = 21, ncol = 5)
for (i in 1:length(genes)) {
    for (j in 1:nrow(ICTYoungMGvsWTMGmatrix)) {
        if (genes[i] == rownames(ICTYoungMGvsWTMGmatrix)[j]) {
            ICTYoungMGvsWTMG[i, ] <- ICTYoungMGvsWTMGmatrix[j, ]
        }
    }
}
rownames(ICTYoungMGvsWTMG) <- relgenes
colnames(ICTYoungMGvsWTMG) <- colnames(ICTYoungMGvsWTMGmatrix)

# Generate gene expression data for Cultured Microglia compared to WT Microglia
topexpressedgenes2 <- topTags(qlf.CultMGvsWTMG,
    n = nrow(qlf.CultMGvsWTMG$table)
)
row_names <- rownames(topexpressedgenes2)
dfrCultMGvsWTMG <- topexpressedgenes2[row_names == genes[1] |
    row_names == genes[2] | row_names == genes[3] | row_names == genes[4] |
    row_names == genes[5] | row_names == genes[6] | row_names == genes[7] |
    row_names == genes[8] | row_names == genes[9] | row_names == genes[10] |
    row_names == genes[11] | row_names == genes[12] | row_names == genes[13] |
    row_names == genes[14] | row_names == genes[15] | row_names == genes[16] |
    row_names == genes[17] | row_names == genes[18] | row_names == genes[19] |
    row_names == genes[20] | row_names == genes[21], ]
CultMGvsWTMGmatrix <- as.matrix(dfrCultMGvsWTMG$table)
CultMGvsWTMG <- matrix(rep(NA, times = 21 * 5), nrow = 21, ncol = 5)
for (i in 1:length(genes)) {
    for (j in 1:nrow(CultMGvsWTMGmatrix)) {
        if (genes[i] == rownames(CultMGvsWTMGmatrix)[j]) {
            CultMGvsWTMG[i, ] <- CultMGvsWTMGmatrix[j, ]
        }
    }
}
rownames(CultMGvsWTMG) <- relgenes
colnames(CultMGvsWTMG) <- colnames(CultMGvsWTMGmatrix)

# Heatmap
library(gplots)

# color palette
library(RColorBrewer)

# Generate matrix of FDR values for genes of interest
logfcmatrix <- cbind(
    ICTMGvsWTMG[, 1], ICTYoungMGvsWTMG[, 1],
    ICTCultMGvsWTMG[, 1], CultMGvsWTMG[, 1]
)
cellnote <- matrix(rep(NA, times = 4 * 21), nrow = 21, ncol = 4)
n <- 0.05 # fdr cut off
for (i in 1:nrow(CultMGvsWTMG)) {
    if (ICTMGvsWTMG[i, 5] < n) {
        cellnote[i, 1] <- "*"
    } else {
        cellnote[i, 1] <- NA
    }
    if (ICTYoungMGvsWTMG[i, 5] < n) {
        cellnote[i, 2] <- "*"
    } else {
        cellnote[i, 2] <- NA
    }
    if (ICTCultMGvsWTMG[i, 5] < n) {
        cellnote[i, 3] <- "*"
    } else {
        cellnote[i, 3] <- NA
    }
    if (CultMGvsWTMG[i, 5] < n) {
        cellnote[i, 4] <- "*"
    } else {
        cellnote[i, 4] <- NA
    }
}

# Heatmap in figure 1C - REPRODUCED ####
pdf(file = "../results/Fig1C.pdf")
heatmap.2(logfcmatrix,
    Colv = NA, Rowv = NA,
    col = rev(brewer.pal(n = 7, name = "RdYlBu")),
    dendrogram = "none", colsep = c(1, 2), rowsep = c(1:21),
    sepcolor = "white", sepwidth = c(0.001, 0.01), trace = "none",
    cexCol = 1, srtCol = 45,
    labCol = c("ICT Adult MG", "ICT P5 MG", "ICT Cultured MG", "Cultured MG"),
    colRow = c(
        rep("Blue", times = 12), rep("Orange", times = 5),
        rep("Red", times = 4)
    ),
    density.info = "none",
    cellnote = cellnote,
    notecex = 2,
    notecol = "black",
    keysize = 1.25, key.xlab = "log2(FC vs WT)", key.title = NA,
    margins = c(12, 9)
)
dev.off()

# volcanoPlot
qlf.ICTCultMGvsWTMG <- glmTreat(fit,
    contrast = c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1), lfc = log2(1.01)
)
qlf.ICTMGvsWTMG <- glmTreat(fit,
    contrast = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1), lfc = log2(1.01)
)
qlf.ICTYoungMGvsWTMG <- glmTreat(fit,
    contrast = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1), lfc = log2(1.01)
)
qlf.CultMGvsWTMG <- glmTreat(fit,
    contrast = c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1), lfc = log2(1.01)
)

# Generate FDR and LogFC data for different conditions
sigICTCultMG <- topTags(qlf.ICTCultMGvsWTMG,
    n = nrow(qlf.ICTCultMGvsWTMG$table)
)
sigICTMG <- topTags(qlf.ICTMGvsWTMG,
    n = nrow(qlf.ICTMGvsWTMG$table)
)
sigICTYoungMG <- topTags(qlf.ICTYoungMGvsWTMG,
    n = nrow(qlf.ICTYoungMGvsWTMG$table)
)
sigCultMG <- topTags(qlf.CultMGvsWTMG, n = nrow(qlf.CultMGvsWTMG$table))

# Plot Volcano Plot of each condition overlapping with Cultured Microglia as in
# Figure 1D - NOT REPRODUCED #####
pdf(file = "../results/Fig1D.pdf", width = 9, height = 3)
par(mfrow = c(1, 3))
with(sigCultMG$table, plot(logFC, -log10(FDR),
    xlab = "log2 (FC vs WT)",
    main = "ICT ADULT MG", pch = 20, col = "orange"
))
points(sigICTMG$table$logFC, -log10(sigICTMG$table$FDR), col = "blue")
with(sigCultMG$table, plot(logFC, -log10(FDR),
    xlab = "log2 (FC vs WT)",
    main = "ICT P5", pch = 20, col = "orange"
))
points(sigICTYoungMG$table$logFC, -log10(sigICTYoungMG$table$FDR),
    col = "blue"
)
with(sigCultMG$table, plot(logFC, -log10(FDR),
    xlab = "log2 (FC vs WT)", main = "ICT CULTURED", pch = 20,
    col = "orange"
))
points(sigICTCultMG$table$logFC, -log10(sigICTCultMG$table$FDR), col = "blue")
dev.off()

#### START OF DATA FOR LAST 3 FIGURES ####

# Combine ICT adult MG, ICT Young MG and ICT Cultured MG into pooled ICT MG
# group
group2 <- as.vector(group1)
group2[c(15, 16, 20, 21, 28, 29, 30, 31, 38, 39, 29)] <- "pICTMG"
group2 <- as.factor(group2)
design2 <- model.matrix(~ 0 + group2, data = RG1$samples)
colnames(design2) <- levels(group2)
RG1 <- estimateDisp(RG1, design2, robust = TRUE)

# Differential Genes Expression
fit2 <- glmQLFit(RG1, design2, robust = TRUE)

qlf.pICTMG <- glmTreat(fit2,
    contrast = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1),
    lfc = log2(1.22)
)
qlf.BMTKO <- glmTreat(fit,
    contrast = c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1), lfc = log2(1.5)
) # BMTKO refers to IP BM in the paper
qlf.ICTBlood <- glmTreat(fit,
    contrast = c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1), lfc = log2(1.5)
)
qlf.ICTBM <- glmTreat(fit,
    contrast = c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1), lfc = log2(1.5)
)
qlf.ICTFetBrain <- glmTreat(fit,
    contrast = c(0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1), lfc = log2(2)
)
qlf.ICTFetLiv <- glmTreat(fit,
    contrast = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1), lfc = log2(1.5)
)
qlf.ICTYS <- glmTreat(fit,
    contrast = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1), lfc = log2(1.5)
)

# Generate gene expression data for pooled ICT Microglia compared to WT
# Microglia
topexpressedgenes <- topTags(qlf.pICTMG, n = nrow(qlf.pICTMG$table))
row_names <- rownames(topexpressedgenes)
dfrpICTMG <- topexpressedgenes[row_names == genes[1] | row_names == genes[2] |
    row_names == genes[3] | row_names == genes[4] | row_names == genes[5] |
    row_names == genes[6] | row_names == genes[7] | row_names == genes[8] |
    row_names == genes[9] | row_names == genes[10] | row_names == genes[11] |
    row_names == genes[12] | row_names == genes[13] | row_names == genes[14] |
    row_names == genes[15] | row_names == genes[16] | row_names == genes[17] |
    row_names == genes[18] | row_names == genes[19] | row_names == genes[20] |
    row_names == genes[21], ]
pICTMGmatrix <- as.matrix(dfrpICTMG$table)
pICTMG <- matrix(rep(NA, times = 21 * 5), nrow = 21, ncol = 5)
for (i in 1:length(genes)) {
    for (j in 1:nrow(pICTMGmatrix)) {
        if (genes[i] == rownames(pICTMGmatrix)[j]) {
            pICTMG[i, ] <- pICTMGmatrix[j, ]
        }
    }
}
rownames(pICTMG) <- relgenes
colnames(pICTMG) <- colnames(pICTMGmatrix)

# Generate gene expression data for IP BM compared to WT Microglia
topexpressedgenes <- topTags(qlf.BMTKO, n = nrow(qlf.BMTKO$table))
row_names <- rownames(topexpressedgenes)
dfrBMTKO <- topexpressedgenes[row_names == genes[1] | row_names == genes[2] |
    row_names == genes[3] | row_names == genes[4] | row_names == genes[5] |
    row_names == genes[6] | row_names == genes[7] | row_names == genes[8] |
    row_names == genes[9] | row_names == genes[10] | row_names == genes[11] |
    row_names == genes[12] | row_names == genes[13] | row_names == genes[14] |
    row_names == genes[15] | row_names == genes[16] | row_names == genes[17] |
    row_names == genes[18] | row_names == genes[19] | row_names == genes[20] |
    row_names == genes[21], ]
BMTKOmatrix <- as.matrix(dfrBMTKO$table)
BMTKO <- matrix(rep(NA, times = 21 * 5), nrow = 21, ncol = 5)
for (i in 1:length(genes)) {
    for (j in 1:nrow(BMTKOmatrix)) {
        if (genes[i] == rownames(BMTKOmatrix)[j]) {
            BMTKO[i, ] <- BMTKOmatrix[j, ]
        }
    }
}
rownames(BMTKO) <- relgenes
colnames(BMTKO) <- colnames(BMTKOmatrix)

# Generate gene expression data for ICT Blood compared to WT Microglia
topexpressedgenes <- topTags(qlf.ICTBlood, n = nrow(qlf.ICTBlood$table))
row_names <- rownames(topexpressedgenes)
dfrICTBlood <- topexpressedgenes[row_names == genes[1] |
    row_names == genes[2] | row_names == genes[3] | row_names == genes[4] |
    row_names == genes[5] | row_names == genes[6] | row_names == genes[7] |
    row_names == genes[8] | row_names == genes[9] | row_names == genes[10] |
    row_names == genes[11] | row_names == genes[12] | row_names == genes[13] |
    row_names == genes[14] | row_names == genes[15] | row_names == genes[16] |
    row_names == genes[17] | row_names == genes[18] | row_names == genes[19] |
    row_names == genes[20] | row_names == genes[21], ]
ICTBloodmatrix <- as.matrix(dfrICTBlood$table)
ICTBlood <- matrix(rep(NA, times = 21 * 5), nrow = 21, ncol = 5)
for (i in 1:length(genes)) {
    for (j in 1:nrow(ICTBloodmatrix)) {
        if (genes[i] == rownames(ICTBloodmatrix)[j]) {
            ICTBlood[i, ] <- ICTBloodmatrix[j, ]
        }
    }
}
rownames(ICTBlood) <- relgenes
colnames(ICTBlood) <- colnames(ICTBloodmatrix)

# Generate gene expression data for ICT BM compared to WT Microglia
topexpressedgenes <- topTags(qlf.ICTBM, n = nrow(qlf.ICTBM$table))
row_names <- rownames(topexpressedgenes)
dfrICTBM <- topexpressedgenes[row_names == genes[1] | row_names == genes[2] |
    row_names == genes[3] | row_names == genes[4] | row_names == genes[5] |
    row_names == genes[6] | row_names == genes[7] | row_names == genes[8] |
    row_names == genes[9] | row_names == genes[10] | row_names == genes[11] |
    row_names == genes[12] | row_names == genes[13] | row_names == genes[14] |
    row_names == genes[15] | row_names == genes[16] | row_names == genes[17] |
    row_names == genes[18] | row_names == genes[19] | row_names == genes[20] |
    row_names == genes[21], ]
ICTBMmatrix <- as.matrix(dfrICTBM$table)
ICTBM <- matrix(rep(NA, times = 21 * 5), nrow = 21, ncol = 5)
for (i in 1:length(genes)) {
    for (j in 1:nrow(ICTBMmatrix)) {
        if (genes[i] == rownames(ICTBMmatrix)[j]) {
            ICTBM[i, ] <- ICTBMmatrix[j, ]
        }
    }
}
rownames(ICTBM) <- relgenes
colnames(ICTBM) <- colnames(ICTBMmatrix)

# Generate gene expression data for ICT Fetal Brain compared to WT Microglia
topexpressedgenes <- topTags(qlf.ICTFetBrain, n = nrow(qlf.ICTFetBrain$table))
row_names <- rownames(topexpressedgenes)
dfrICTFetBrain <- topexpressedgenes[row_names == genes[1] |
    row_names == genes[2] | row_names == genes[3] | row_names == genes[4] |
    row_names == genes[5] | row_names == genes[6] | row_names == genes[7] |
    row_names == genes[8] | row_names == genes[9] | row_names == genes[10] |
    row_names == genes[11] | row_names == genes[12] | row_names == genes[13] |
    row_names == genes[14] | row_names == genes[15] | row_names == genes[16] |
    row_names == genes[17] | row_names == genes[18] | row_names == genes[19] |
    row_names == genes[20] | row_names == genes[21], ]
ICTFetBrainmatrix <- as.matrix(dfrICTFetBrain$table)
ICTFetBrain <- matrix(rep(NA, times = 21 * 5), nrow = 21, ncol = 5)
for (i in 1:length(genes)) {
    for (j in 1:nrow(ICTFetBrainmatrix)) {
        if (genes[i] == rownames(ICTFetBrainmatrix)[j]) {
            ICTFetBrain[i, ] <- ICTFetBrainmatrix[j, ]
        }
    }
}
rownames(ICTFetBrain) <- relgenes
colnames(ICTFetBrain) <- colnames(ICTFetBrainmatrix)

# Generate gene expression data for ICT Fetal Liver compared to WT Microglia
topexpressedgenes <- topTags(qlf.ICTFetLiv, n = nrow(qlf.ICTFetLiv$table))
row_names <- rownames(topexpressedgenes)
dfrICTFetLiv <- topexpressedgenes[row_names == genes[1] |
    row_names == genes[2] | row_names == genes[3] | row_names == genes[4] |
    row_names == genes[5] | row_names == genes[6] | row_names == genes[7] |
    row_names == genes[8] | row_names == genes[9] | row_names == genes[10] |
    row_names == genes[11] | row_names == genes[12] | row_names == genes[13] |
    row_names == genes[14] | row_names == genes[15] | row_names == genes[16] |
    row_names == genes[17] | row_names == genes[18] | row_names == genes[19] |
    row_names == genes[20] | row_names == genes[21], ]
ICTFetLivmatrix <- as.matrix(dfrICTFetLiv$table)
ICTFetLiv <- matrix(rep(NA, times = 21 * 5), nrow = 21, ncol = 5)
for (i in 1:length(genes)) {
    for (j in 1:nrow(ICTFetLivmatrix)) {
        if (genes[i] == rownames(ICTFetLivmatrix)[j]) {
            ICTFetLiv[i, ] <- ICTFetLivmatrix[j, ]
        }
    }
}
rownames(ICTFetLiv) <- relgenes
colnames(ICTFetLiv) <- colnames(ICTFetLivmatrix)

# Generate gene expression data for ICT Yolk sacs compared to WT Microglia
topexpressedgenes <- topTags(qlf.ICTYS, n = nrow(qlf.ICTYS$table))
row_names <- rownames(topexpressedgenes)
dfrICTYS <- topexpressedgenes[row_names == genes[1] | row_names == genes[2] |
    row_names == genes[3] | row_names == genes[4] | row_names == genes[5] |
    row_names == genes[6] | row_names == genes[7] | row_names == genes[8] |
    row_names == genes[9] | row_names == genes[10] | row_names == genes[11] |
    row_names == genes[12] | row_names == genes[13] | row_names == genes[14] |
    row_names == genes[15] | row_names == genes[16] | row_names == genes[17] |
    row_names == genes[18] | row_names == genes[19] | row_names == genes[20] |
    row_names == genes[21], ]
ICTYSmatrix <- as.matrix(dfrICTYS$table)
ICTYS <- matrix(rep(NA, times = 21 * 5), nrow = 21, ncol = 5)
for (i in 1:length(genes)) {
    for (j in 1:nrow(ICTYSmatrix)) {
        if (genes[i] == rownames(ICTYSmatrix)[j]) {
            ICTYS[i, ] <- ICTYSmatrix[j, ]
        }
    }
}
rownames(ICTYS) <- relgenes
colnames(ICTYS) <- colnames(ICTYSmatrix)

# Generate matrix of significant (*) based on FDR values
logfcmatrix1 <- cbind(
    pICTMG[, 1], ICTFetBrain[, 1], ICTYS[, 1],
    ICTFetLiv[, 1], ICTBM[, 1], ICTBlood[, 1], BMTKO[, 1]
)
cellnote <- matrix(rep(NA, times = 7 * 21), nrow = 21, ncol = 7)
n <- 0.05 # fdr cut off
for (i in 1:nrow(ICTBM)) {
    if (pICTMG[i, 5] < n) {
        cellnote[i, 1] <- "*"
    } else {
        cellnote[i, 1] <- NA
    }
    if (ICTFetBrain[i, 5] < n) {
        cellnote[i, 2] <- "*"
    } else {
        cellnote[i, 2] <- NA
    }
    if (ICTYS[i, 5] < n) {
        cellnote[i, 3] <- "*"
    } else {
        cellnote[i, 3] <- NA
    }
    if (ICTFetLiv[i, 5] < n) {
        cellnote[i, 4] <- "*"
    } else {
        cellnote[i, 4] <- NA
    }
    if (ICTBM[i, 5] < n) {
        cellnote[i, 5] <- "*"
    } else {
        cellnote[i, 5] <- NA
    }
    if (ICTBlood[i, 5] < n) {
        cellnote[i, 6] <- "*"
    } else {
        cellnote[i, 6] <- NA
    }
    if (BMTKO[i, 5] < n) {
        cellnote[i, 7] <- "*"
    } else {
        cellnote[i, 7] <- NA
    }
}

# Heatmap for figure 4A - REPRODUCED ####
pdf(file = "../results/Fig4A.pdf")
heatmap.2(logfcmatrix1,
    Colv = NA, Rowv = NA, col = rev(brewer.pal(n = 9, name = "RdYlBu")),
    dendrogram = "none", colsep = c(1:7), rowsep = c(1:21),
    sepcolor = "white", sepwidth = c(0.001, 0.01), trace = "none",
    cexCol = 1, srtCol = 45,
    labCol = c(
        "pICT MG", "ICT Fet Brain", "ICT YS", "ICT Fet Liver",
        "ICT BM", "ICT Blood", "BMTKO"
    ),
    colCol = c(rep("Blue", times = 3), rep("Orange", times = 4)),
    colRow = c(
        rep("Blue", times = 12), rep("Orange", times = 5),
        rep("Red", times = 4)
    ),
    density.info = "none", cellnote = cellnote, notecex = 2,
    notecol = "black", keysize = 1.25,
    key.xlab = "log2 (ICT vs WT)", key.title = NA
)
dev.off()

# Principal Component Analysis Plot
# Import FPKM data
FPKMdata <- read.csv("FPKMall Data.csv",
    header = TRUE, sep = ",",
    quote = "\"'", stringsAsFactors = FALSE
)
rownames(FPKMdata) <- FPKMdata[, 1]
FPKMdata <- FPKMdata[, -1]
FPKMdata_log <- log2(FPKMdata + 1) # find log2(FPKM+1) value for each gene
FPKMdata_logMat <- as.matrix(FPKMdata_log)
gene_variance <- apply(X = FPKMdata_logMat, MARGIN = 1, FUN = function(x) {
    var(x, na.rm = TRUE)
})
FPKMdata$variance <- gene_variance
FPKMord <- FPKMdata[order(-FPKMdata$variance), ] # order the genes based on their variance
FPKMfilt <- as.matrix(FPKMord[1:2500, -53]) # select the top 2500 genes
colnames(FPKMfilt) <- c(
    rep("WTMG", times = 7),
    rep("pICTMG", times = 2),
    rep("pICTMG", times = 4),
    rep("pICTMG", times = 4),
    rep("ICTFetBrain", times = 5),
    rep("ICTYS", times = 3),
    rep("ICTFetLiv", times = 7),
    rep("ICTBM", times = 3),
    rep("ICTBlood", times = 3), rep("Acute11B", times = 3),
    rep("Cult11B", times = 3),
    rep("MBTCtl", times = 3), rep("BMTKO", times = 5)
)
FPKMfilt <- FPKMfilt[, -c(39:41, 45:47)]
FPKMfinal <- t(FPKMfilt) # transpose the matrix so that the conditions are in
# the column and the genes are in the row

# Calculate PCA with SVD imputation and unit variance scaling
### LOAD BCV PACKAGE ###
library(bcv)

SVDimputation <- impute.svd(FPKMfinal, maxite = 10000)$x # using SVD imputation
# to generate missing data
FPKM.pca1 <- prcomp(SVDimputation, center = TRUE, scale. = TRUE) # scale is
# TRUE means that the pca is calculated with unit variance scaling

# Plot PCA component 1 and 2 using ggbiplot - NOT REPRODUCED ####
library(ggbiplot)
groupFPKM <- c(
    rep("WT MG", times = 7), rep("pICT MG", times = 10),
    rep("ICT Fetal Brain", times = 5), rep("ICT YS", times = 3),
    rep("ICT Fetal Liver", times = 7), rep("ICT Blood/BM", times = 6),
    rep("Cultured MG", times = 3), rep("IP BM", times = 5)
)

pdf("../results/Fig4B.pdf")
ggbiplot(FPKM.pca1,
    ellipse = TRUE, ellipse.prob = 0.95, obs.scale = 1,
    groups = groupFPKM, var.axes = FALSE
) +
    scale_colour_manual(
        name = "Condition",
        values = c(
            "Orange", "Violet", "Red", "Blue", "Yellow", "lightblue",
            "Green", "Black"
        )
    ) +
    theme(legend.position = "bottom") +
    scale_y_reverse() +
    scale_x_reverse()
dev.off()

# Generate the 1000 most variant genes based on log2(FPKMaverage +1) value
FPKMavg <- read.csv("FPKMaverage.csv",
    header = TRUE, sep = ",",
    quote = "\"'", stringsAsFactors = FALSE
)
rownames(FPKMavg) <- FPKMavg[, 1]
FPKMavg <- FPKMavg[, -1]
FPKMavg_log <- log2(FPKMavg + 1) # find log2(FPKM+1) value for each gene
FPKMavg_logMat <- as.matrix(FPKMavg_log)
gene_avg_variance <- apply(X = FPKMavg_logMat, MARGIN = 1, FUN = function(x) {
    var(x, na.rm = TRUE)
})
FPKMavg$variance <- gene_avg_variance
FPKMavgord <- FPKMavg[order(-FPKMavg$variance), ]
FPKMavg1000 <- as.matrix(FPKMavgord[1:1000, -c(2:4, 10, 12, 15:17)])
colnames(FPKMavg1000) <- c(
    "WT MG", "ICT Fetal Br", "ICT YS",
    "ICT Fetal Liv", "ICT BM", "ICT Blood", "Cultured MG", "IP BM",
    "pICT MG"
)

# Generate Unsupervised hierarchical clustering
library(pvclust)

# Create a function to calculate corvariance using Spearman method
spearman <- function(x, ...) {
    x <- as.matrix(x)
    res <- as.dist(1 - cor(x, method = "spearman", use = "everything"))
    res <- as.dist(res)
    attr(res, "method") <- "spearman"
    return(res)
}
result <- pvclust(FPKMavg1000,
    method.hclust = "single",
    method.dist = spearman, use.cor = "all.obs",
    nboot = 10000
)

pdf("../results/Fig4C.pdf")
plot(result,
    print.pv = TRUE, print.num = FALSE, float = 0.03,
    col.pv = c("grey", 2)
)
dev.off()
