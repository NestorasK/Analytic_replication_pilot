# Figure 2F - REPRODUCED ####
rm(list = ls())
library("ggplot2")
setwd("PATH_TO_Analytic_replication_pilot")

mt_poly_a <- data.frame(
    tissue = factor(
        x = c(
            "heart", "colon", "muscle", "kidney",
            "adipose", "breast", "liver", "brain", "testes", "thyroid",
            "adrenal", "prostate", "lymph", "ovary", "whiteblood", "lung"
        ),
        levels = c(
            "heart", "colon", "muscle", "kidney", "adipose", "breast", "liver",
            "brain", "testes", "thyroid", "adrenal", "prostate", "lymph",
            "ovary", "whiteblood", "lung"
        )
    ),
    proportion = c(
        0.29, 0.225, 0.195, 0.19, 0.15, 0.145, 0.14, 0.13, 0.10,
        0.06, 0.059, 0.055, 0.054, 0.053, 0.053, 0.052
    )
)
ggplot(data = mt_poly_a, aes(x = tissue, y = proportion)) +
    geom_bar(stat = "identity", color = "black", fill = "lightgrey") +
    labs(
        title = "", x = "",
        y = "mitochondrial contribution\nto cellular polyA+ RNA",
        ylim = c(0, 3.5)
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(
    filename = "team1/results/Figure_2F.png",
    width = 4, height = 4
)
