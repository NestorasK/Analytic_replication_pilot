rm(list = ls())
library(ggplot2)
setwd("PATH_TO_Analytic_replication_pilot")

## Generate a csv table first and read in the R console-generates a dataframe #

table3 <- read.csv(
    file = "team3/data/Fig3d_e_f.csv", header = TRUE,
    stringsAsFactors = TRUE, skip = 1
)

### remove rows not defined by samples/dicer_mutation
samples <- c(
    "TCGA_CHOL", "TCGA_COAD", "TCGA_UCEC", "TCGA_SKCM",
    "TCGA_THCA", "TCGA_UCS"
)

mutation_types <- c(
    "R944Q", "RNAseIIIB_hotspot_biallelic",
    "RNAseIIIB_hotspot", "RNAseIIIB_S1344L_biallelic"
)

table3new <- table3[table3$TCGA_NAME %in% samples, ]

### Number of samples by cancer_name_Here when I  run this code the number of
# samples for each mutation type in figure 3f matches for some cancer while
# not for others??
num_TCGA_NAME <- data.frame(table(table3new$TCGA_NAME))
colnames(num_TCGA_NAME)[1] <- "TCGA_NAME"


table3fnew <- merge(
    x = table3new, y = num_TCGA_NAME,
    by = "TCGA_NAME"
)
table3fnew$TCGA_NAME2Plot <- paste(
    (table3fnew$TCGA_NAME), "(",
    table3fnew$Freq, ")"
)

boxplotTCGA <- as.data.frame(table3fnew[, c(7, 3)])

### Plot based on TCGA_names

plot_i <- qplot(
    x = table3fnew$TCGA_NAME2Plot,
    y = table3fnew$log2.miR.5p.miR.3p.,
    col = factor(droplevels(
        x = table3fnew$Mutation_type,
        exclude = c("WT", "DICER_other"), pch = 1
    )),
    xlab = ""
) +
    theme_classic() +
    geom_jitter(width = 0.3) +
    geom_hline(yintercept = 0, col = "red", lwd = 0.5, lty = 2) +
    ylab("log2(miR-5p/miR-3p)") +
    ylim(-1.5, 3) +
    theme(
        axis.title = element_text(size = 10),
        axis.text.x = element_text(angle = 60, hjust = 1)
    ) +
    theme(legend.title = element_blank()) +
    geom_boxplot(
        data = boxplotTCGA,
        aes(
            x = factor(boxplotTCGA$TCGA_NAME),
            y = boxplotTCGA$log2.miR.5p.miR.3p.
        ), alpha = 0.01,
        inherit.aes = FALSE, col = "black"
    )
ggsave(
    filename = "team3/results/Figure_3F.pdf", plot = plot_i,
    width = 7, height = 5
)
