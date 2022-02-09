#### 3E###
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

#### number of samples by mutation type (??Why does ovarian tumor still show up)
num_mut_type <- data.frame(table(table3new$Mutation_type))
num_mut_type

colnames(num_mut_type)[1] <- "Mutation_type"

table3e <- merge(
      x = table3new, y = num_mut_type,
      by = "Mutation_type"
)
table3e$Mutation_type2Plot <- paste(
      (table3e$Mutation_type),
      "(", table3e$Freq, ")"
)

## qplot and save in pdf(commented out)
plot_i <- qplot(
      x = table3e$Mutation_type2Plot,
      y = table3e$log2.miR.5p.miR.3p.,
      col = factor(x = table3e$Mutation_type2Plot),
      xlab = ""
) +
      theme_classic() +
      ylab("log2(miR-5p/miR-3p)") +
      ylim(-1.5, 3) +
      theme(
            axis.title = element_text(size = 10),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            axis.ticks.x = element_blank()
      ) +
      theme(legend.title = element_blank()) +
      geom_jitter(width = 0.2) +
      geom_hline(yintercept = 0, col = "red", lwd = 0.5, lty = 2) +
      geom_boxplot(
            color = c(rep("gray", times = 5), "black"),
            alpha = 0.01
      )

ggsave(
      filename = "team3/results/Figure_3E.pdf", plot = plot_i,
      width = 7, height = 5
)
