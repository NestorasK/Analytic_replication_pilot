#### Fig3D###
rm(list = ls())
library(ggplot2)
setwd("PATH_TO_Analytic_replication_pilot")

## Load the csv file onto R console with approapriate
data <- read.csv(file = "team3/data/Fig3d_e_f.csv", header = TRUE, skip = 1)
data

fig3d <- qplot(
    y = sort(data$log2.miR.5p.miR.3p., decreasing = TRUE), ylim = c(-1, 3),
    ylab = "mi53=log2(miR-5p/miR-3p)"
) +
    theme_classic() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    xlab("9919 TCGA samples") +
    geom_area()

ggsave(
    filename = "team3/results/Figure_3D.pdf", plot = fig3d,
    width = 4, height = 3
)
