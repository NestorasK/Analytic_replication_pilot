# Figure 2B - NOT REPRODUCED #####
rm(list = ls())
library("ggplot2")
setwd("PATH_TO_Analytic_replication_pilot")

rna_transcripts <-
    data.frame(
        RNA = c(
            "rRNA", "tRNA", "mRNA",
            "intergenic/antisense"
        ),
        proportion = c(78, 13, 8, 1)
    )
ggplot(rna_transcripts, aes(x = "", y = proportion, fill = RNA)) +
    geom_bar(width = 1, stat = "identity", color = "black") +
    coord_polar("y", start = 0, direction = -1) +
    labs(title = "", x = "", y = "") +
    scale_fill_manual(values = c(
        "grey", "royalblue2",
        "indianred2", "seagreen2"
    )) +
    theme_void() +
    geom_text(aes(y = proportion, label = proportion),
        color = "black",
        size = 5, position = position_stack(0.5)
    )
ggsave(
    filename = "team1/results/Figure_2B.png",
    width = 4,
    height = 4
)
