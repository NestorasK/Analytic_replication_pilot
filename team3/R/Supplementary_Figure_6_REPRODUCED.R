rm(list = ls())
library(ggplot2)
setwd("PATH_TO_Analytic_replication_pilot")

## Generate a csv table first and read in the R console-generates a dataframe##
table3 <- read.csv(
      file = "team3/data/Fig3d_e_f.csv", header = TRUE,
      stringsAsFactors = TRUE, skip = 1
)

### remove rows not defined by samples/dicer_mutation
samples <- c(
      "WT", "R944Q", "DICER_other", "RNAseIIIB_hotspot_biallelic",
      "RNAseIIIB_hotspot", "RNAseIIIB_S1344L_biallelic"
)
table3new <- table3[table3$Mutation_type %in% samples, ]


### Number of samples by cancer_name
num_TCGA_NAME <- data.frame(table(table3new$TCGA_NAME))
colnames(num_TCGA_NAME)[1] <- "TCGA_NAME"

table3f <- merge(
      x = table3new, y = num_TCGA_NAME,
      by = "TCGA_NAME"
)
table3f$TCGA_NAME2Plot <- paste((table3f$TCGA_NAME), "(", table3f$Freq, ")")


boxplotTCGA <- as.data.frame(table3f[, c(1, 3)])



### Plot based on TCGA_names
plot_i <- qplot(
      x = table3f$TCGA_NAME,
      y = table3f$log2.miR.5p.miR.3p.,
      col = factor(droplevels(
            x = table3f$Mutation_type,
            exclude = c("WT", "DICER_other"), pch = 1
      )),
      xlab = ""
) +
      theme_classic() +
      ylab("log2(miR-5p/miR-3p)") +
      ylim(-1.5, 3) +
      theme(
            axis.title = element_text(size = 10),
            axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      theme(legend.title = element_blank()) +
      geom_jitter(width = 0.15) +
      geom_hline(yintercept = 0, col = "red", lwd = 0.5, lty = 2) +
      geom_boxplot(data = boxplotTCGA, aes(
            x = factor(boxplotTCGA$TCGA_NAME),
            y = boxplotTCGA$log2.miR.5p.miR.3p.
      ), alpha = 0.01, inherit.aes = FALSE, col = "black")

ggsave(
      filename = "team3/results/Supplementary_Figure6.pdf",
      plot = plot_i, width = 10, height = 10
)
