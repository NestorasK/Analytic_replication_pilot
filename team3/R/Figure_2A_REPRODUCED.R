##### Figure 2 Cancer DICER PAPER####
# set your working directory to team3 folder
# setwd("PATH_TO/team3")
rm(list = ls())
library(ggplot2)
library(ggpubr)
library(reshape2)

setwd("PATH_TO_Analytic_replication_pilot")

### Read the .csv file from correct directory to the global environment
tcga <- read.table(
    file = "team3/data/tcga_clinical.csv",
    header = TRUE, sep = ",", stringsAsFactors = FALSE
)

### Count Total Samples by Cancer_TYPE_ACRONYM
sample_count <- table(tcga$CANCER_TYPE_ACRONYM)

###### Count Total Samples by Cancer_TYPE_ACRONYM which excludes none factor
CANCER_DICER <- table(tcga$CANCER_TYPE_ACRONYM[tcga$dicer_muts != "none"])

########## Count total samples which have no mutations by Cancer Type
No_CANCER <- table(tcga$CANCER_TYPE_ACRONYM[tcga$Dicer_Mut_Cnt == 0])

####### Count total samples which has some mutations by Cancer Type
total_Cancer_TYPE <- sample_count - No_CANCER

## Calculate the ratio of samples with DICER 1 Mutations by CANCER SubTYPE
ratio <- ((total_Cancer_TYPE) / sample_count) * 100

### Create an table that has os in place for all kinds of dicer mutations
table_dicerfreq <- matrix(
    data = 0,
    nrow = length(unique(tcga$CANCER_TYPE_ACRONYM)),
    ncol = length(unique(tcga$dicer_muts))
)

### Replace the table row  names and col names by CANCER TYPE and types of
# DICER mutations respectively
rownames(table_dicerfreq) <- unique(tcga$CANCER_TYPE_ACRONYM)
colnames(table_dicerfreq) <- unique(tcga$dicer_muts)


### Generate the for looop to count dicer_mutaions by unique CANCER_TYPE
for (i in 1:length(rownames(table_dicerfreq))) {
    cancer_i <- rownames(table_dicerfreq)[i]
    CANCER_TYPE_i <- tcga[tcga$CANCER_TYPE_ACRONYM == cancer_i, ]
    CANCER_TYPE_i
    dicer_muts_i <- table(CANCER_TYPE_i$dicer_muts)
    dicer_muts_i

    #### Frequency of dicer_mutations by 1st CANCER_SYBTYPE
    dicer_muts_i.fre <- (dicer_muts_i /
        total_Cancer_TYPE[names(total_Cancer_TYPE) == cancer_i])
    dicer_muts_i.fre

    #### Generate an index for columns where names of elements in
    # vector(dicer_muts_i.fre) matches colnames in table_dicerfreq
    ind <- match(names(dicer_muts_i.fre), colnames(table_dicerfreq))
    ind
    table_dicerfreq[i, ind] <- dicer_muts_i.fre
    table_dicerfreq[i, ind]
}
table_dicerfreq <- table_dicerfreq[, -2] ## remove the none column

### Consolidate ratio table and table_dicer_freq -possible because rows match ,
# otherwise figure out how to use merge function to combine two tables
# that have same rows in unmatched order

merged_table <- cbind(table_dicerfreq, ratio)

### Coerce the data into a dataframe
merged_table <- as.data.frame(merged_table)

# Normalize the dicermutationtype for each cancer's total mutation frequency
merged_table[1:5] <- lapply(merged_table[1:5], "*", merged_table$ratio)


merged_table$CancerType <- row.names(merged_table)
row.names(merged_table) <- NULL

### Melt the table by ratio/yvalue and cancer_type_acronynm

melted_table <- melt(
    data = merged_table,
    id.vars = c("ratio", "CancerType"),
    variable.name = "MutationType",
    value.name = "Height"
)

### Add a column to the melted table with number of samples (recycles the cancer type)

melted_table$CancerType_numSample <- paste(names(sample_count),
    "(n=", sample_count, ")",
    sep = ""
)


### Generate an index by decreasing ratio from the melted table
index <- order(melted_table$ratio, decreasing = TRUE)


#### Plot the data using index to sort x-axis(CANCERTYPE) based on
# sorted y-axis(% of total mutations)


plot_i <- qplot(
    x = factor(melted_table$CancerType_numSample[index],
        levels = unique(melted_table$CancerType_numSample[index])
    ),
    y = melted_table$Height[index],
    ylab = "",
    xlab = "",
    geom = "blank",
    fill = melted_table$MutationType[index]
) +
    theme_classic() +
    geom_bar(stat = "identity") +
    ylab("% Dicer Mutation(per tumor subtype)") +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(legend.position = c(0.7, 0.7), legend.title = element_blank()) +
    scale_y_continuous(
        breaks = c(0, 5, 10, 25, 50, 75, 100),
        limits = c(0, 75)
    )
ggsave(
    filename = "team3/results/Figure_2A.pdf", plot = plot_i,
    width = 12, height = 8
)
