# Figure 2A - NOT REPRODUCED #####
rm(list = ls())
library("ggplot2")
setwd("PATH_TO_Analytic_replication_pilot")

promotertranscripts <- data.frame(promoter = c(
    "HSP (major)",
    "HSP (minor)", "LSP"
), proportion = c(80, 9, 11))
ggplot(promotertranscripts, aes(x = "", y = proportion, fill = promoter)) +
    geom_bar(width = 1, stat = "identity", color = "black") +
    coord_polar("y", start = 0, direction = -1) +
    labs(title = "", x = "", y = "") +
    scale_fill_manual(values = c("steelblue1", "royalblue2", "grey")) +
    theme_void() +
    geom_text(aes(y = proportion, label = proportion),
        color = "black",
        size = 5, position = position_stack(0.5)
    )
ggsave(
    filename = "team1/results/Figure_2A.png",
    height = 3,
    width = 3
)
