# Figure 5A - NOT REPRODUCED #####
rm(list = ls())
library("ggplot2")
setwd("PATH_TO_Analytic_replication_pilot")

# Circular and pie chart in R
library(circlize)
spread <- read.csv(
    file = "team1/data/Figure5A.csv", header = TRUE,
    stringsAsFactors = TRUE
)
mt_genome <- read.csv(
    file = "team1/data/MtGenome.csv", header = TRUE,
    stringsAsFactors = TRUE
)
mt_genome_h <- mt_genome[mt_genome$Strand == "+", c(1, 2, 4)]
mt_genome_l <- mt_genome[mt_genome$Strand == "-", c(1, 2, 4)]

# For the heavy strand (+)
heavy_strand <- spread[spread$Heavy......Light.... == "+", c(1, 2, 3, 5)]
circos.genomicInitialize(heavy_strand, plotType = NULL)
circos.trackPlotRegion(
    ylim = c(0, 17000),
    bg.col = "lightblue",
    bg.border = NA, track.height = 0.05
)
circos.genomicLines(
    region = heavy_strand[, c(2, 3)], value = heavy_strand[, c(2, 3)],
    type = "h", col = "red"
)
circos.text(
    x = mt_genome_h[, 1], y = mt_genome_h[, 2], labels = mt_genome_h[, 3],
    niceFacing = TRUE,
    adj = c(1, 0.5), facing = "reverse.clockwise", cex = 0.5
)

# For the light strand (-)
light_strand <- spread[spread$Heavy......Light.... == "-", c(1, 2, 3, 5)]
circos.genomicInitialize(light_strand, plotType = NULL)
circos.trackPlotRegion(
    ylim = c(0, 17000),
    bg.col = "lightgrey", bg.border = "lightgrey", track.height = 0.05
)
circos.genomicLines(
    region = light_strand[, c(2, 3)], value = light_strand[, c(2, 3)],
    type = "h", col = "red"
)
circos.text(
    x = mt_genome_l[, 1], y = mt_genome_l[, 2],
    labels = mt_genome_l[, 3], niceFacing = TRUE,
    adj = c(1, 0.5), facing = "clockwise", cex = 0.5,
)


# To put on the same plot
pdf(
    file = "team1/results/Figure5A_circos.pdf",
    width = 7, height = 7
)
circos.genomicInitialize(heavy_strand, plotType = NULL)
circos.trackPlotRegion(
    ylim = c(0, 17000),
    bg.col = "lightblue",
    bg.border = NA, track.height = 0.05
)
circos.genomicLines(
    region = heavy_strand[, c(2, 3)], value = heavy_strand[, c(2, 3)],
    type = "h", col = "red"
)
circos.text(
    x = mt_genome_h[, 1], y = mt_genome_h[, 2],
    labels = mt_genome_h[, 3], niceFacing = TRUE,
    adj = c(1, 0.5), facing = "reverse.clockwise", cex = 0.5
)

circos.trackPlotRegion(
    ylim = c(0, 17000),
    bg.col = "lightgrey", bg.border = "lightgrey", track.height = 0.05
)
circos.genomicLines(
    region = light_strand[, c(2, 3)], value = light_strand[, c(2, 3)],
    type = "h", col = "red"
)
circos.text(
    x = mt_genome_l[, 1], y = mt_genome_l[, 2],
    labels = mt_genome_l[, 3], niceFacing = TRUE,
    adj = c(1, 0.5), facing = "clockwise", cex = 0.5,
)
dev.off()

# For the pie chart
small_rna <- data.frame(
    RNA = c("tRNAs", "mRNAs", "rRNAs"),
    proportion = c(64, 28, 8)
)
ggplot(small_rna, aes(x = "", y = proportion, fill = RNA)) +
    geom_bar(width = 1, stat = "identity", color = "black") +
    coord_polar("y", start = 0, direction = -1) +
    labs(title = "proportion of\nsmall RNAs", x = "", y = "") +
    scale_fill_manual(values = c("lightblue", "green", "darkblue")) +
    theme_void() +
    geom_text(aes(y = proportion, label = proportion),
        color = "white", size = 5, position = position_stack(0.5)
    )
ggsave(
    filename = "team1/results/Figure5A_pie_chart.png",
    width = 3, height = 3
)
