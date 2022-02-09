# Figure 2C - REPRODUCED ####
# Data available in file: 1-s2.0-S0092867411007677-mmc2.xls, sheet A
# of the initial publication
rm(list = ls())
library("ggplot2")
library("gplots")
setwd("PATH_TO_Analytic_replication_pilot")

mitoplastexpression <-
    data.frame(
        RNA = c(
            "12S", "12S", "16S", "16S", "ATPase8/6",
            "ATPase8/6", "CO1", "CO1", "CO2", "CO2", "CO3", "CO3",
            "Cytb", "Cytb", "ND1", "ND1", "ND2", "ND2", "ND3", "ND3",
            "ND4L/4", "ND4L/4", "ND5", "ND5", "ND6", "ND6"
        ),
        Strand = rep(x = c("sense", "antisense"), times = 13, each = 1),
        RPKM = c(
            22.3566, 0.0138, 39.5711, 0.0051, 0.5768, 0.0375, 1.5706,
            0.0253, 1.6935, 0.0315, 0.5922, 0.0485, 0.3328, 0.1443, 0.4694,
            0.0340, 0.6519, 0.0214, 0.8257, 0.0973, 0.5409, 0.0723, 0.2293,
            0.2095, 0.1272, 0.0958
        ),
        CI = c(
            0.0907, 0.0007, 0.1207, 0.0004, 0.0107, 0.001, 0.0240, 0.0009,
            0.025, 0.001, 0.0106, 0.0012, 0.0111, 0.0021, 0.0131, 0.001,
            0.0155, 0.0008, 0.0174, 0.0018, 0.0141, 0.0015, 0.0092,
            0.0026, 0.002, 0.0059
        )
    )

# First plot with y axis showing values 0-40 rpkm^104 (12S and 16S expression)
ggplot(data = mitoplastexpression, aes(x = RNA, y = RPKM, fill = Strand)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = "", x = "", y = "expression (rpkm x 10^4)") +
    scale_fill_manual(values = c("darkgreen", "blue")) +
    geom_errorbar(aes(ymin = RPKM - CI, ymax = RPKM + CI),
        width = .2, position = position_dodge(0.9)
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(
    filename = "team1/results/Figure_2C_1.png",
    width = 4, height = 4
)

# Second plot with y axis showing values 0-2 rpkm^104
# (Not including 12S and 16S)
mitoplastexpression2 <- mitoplastexpression[-c(1:4), ]
ggplot(data = mitoplastexpression2, aes(x = RNA, y = RPKM, fill = Strand)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = "", x = "", y = "expression (rpkm x 10^4)", ylim = c(0, 2)) +
    scale_fill_manual(values = c("darkgreen", "blue")) +
    geom_errorbar(aes(ymin = RPKM - CI, ymax = RPKM + CI),
        width = .2,
        position = position_dodge(0.9)
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = "team1/results/Figure_2C_2.png", width = 4, height = 4)
