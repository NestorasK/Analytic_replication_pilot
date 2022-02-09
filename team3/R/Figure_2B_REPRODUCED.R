##### Figure 2 Cancer DICER PAPER####
rm(list = ls())
library(ggplot2)
library(ggpubr)
library(reshape2)

setwd("PATH_TO_Analytic_replication_pilot")


### Read the .csv file from correct directory to the global environment

msk <- read.csv(
    file = "team3/data/TableS2_IMPACT_DICER1_mutants_4252019.csv",
    header = TRUE, sep = ",", stringsAsFactors = FALSE,
    skip = 1
)

### Count Total Samples by Cancer_TYPE_ACRONYM retaining samples with more
# than 50 cases
sample_count <- table(msk$CANCER_TYPE)

### Cancer types which has more than 50 cases
sample_50 <- sample_count[sample_count > 50]

## Extract names for cancers which has more than 50 cases
CANCER_50 <- names(x = sample_50)


## Add the two hyper_mutated /POLE mutation cancers which are present in
# TCGA data sets.

# CANCER_all<-c(CANCER_50,c("Uterine Mixed Endometrial Carcinoma",
#  "Colorectal Cancer Hypermutated"))
## when we do this there are still more cancer_types that fit the criteria
# (61-mine for, 41 from their description in the paper)
### They take various sub-types of cancer and coalesce for some types while not
# for others and I could not find a clear methods to be able to recreate the
# figure. Went ahead anyways.


### Subset only cancers that have more than 50 cases ,
# 59 cancer types without the POLE mutations included
msk_50 <- msk[msk$CANCER_TYPE %in% CANCER_50, ]

###
sample_count50 <- table(msk_50$CANCER_TYPE)


###### Count Total Samples by Cancer_TYPE_ACRONYM which excludes none factor
CANCER_DICER <- table(msk_50$CANCER_TYPE[msk_50$dicer_muts != "none"])



########## Count total samples which have no mutations by Cancer Type
No_CANCER <- table(msk_50$CANCER_TYPE[msk_50$number_of_mutations == 0])


####### Count total samples which has some mutations by Cancer Type
total_Cancer_TYPE <- sample_count50 - No_CANCER


## Calculate the ratio of samples with DICER 1 Mutations by CANCER SubTYPE
ratio <- ((total_Cancer_TYPE) / sample_count50) * 100


### Create an table that has os in place for all kinds of dicer mutations

table_dicerfreq <- matrix(
    data = 0,
    nrow = length(unique(msk_50$CANCER_TYPE)),
    ncol = length(unique(msk_50$dicer_muts))
)


### Replace the table row  names and col names by CANCER TYPE and types of DICER
# mutations respectively

rownames(table_dicerfreq) <- unique(msk_50$CANCER_TYPE)
colnames(table_dicerfreq) <- unique(msk_50$dicer_muts)


### Generate the for looop to count dicer_mutaions by unique CANCER_TYPE

for (i in 1:length(rownames(table_dicerfreq))) {
    cancer_i <- rownames(table_dicerfreq)[i]
    CANCER_TYPE_i <- msk_50[msk_50$CANCER_TYPE == cancer_i, ]
    CANCER_TYPE_i
    dicer_muts_i <- table(CANCER_TYPE_i$dicer_muts)
    dicer_muts_i

    #### Frequency of dicer_mutations by 1st CANCER_SUBTYPE
    dicer_muts_i.fre <- (dicer_muts_i / total_Cancer_TYPE[
        names(total_Cancer_TYPE) == cancer_i
    ])
    dicer_muts_i.fre

    #### Generate an index for columns where names of elements in vector
    # (dicer_muts_i.fre) matches colnames in table_dicerfreq
    ind <- match(names(dicer_muts_i.fre), colnames(table_dicerfreq))
    ind
    table_dicerfreq[i, ind] <- dicer_muts_i.fre
    table_dicerfreq[i, ind]
}
table_dicerfreq <- table_dicerfreq[, -2] ## remove the none column



### Consolidate ratio table and table_dicer_freq -possible because rows match ,
# otherwise figure out how to use merge function to combine two tables that have
# same rows in unmatched order

merged_table <- cbind(table_dicerfreq, ratio)

### Coerce the data into a dataframe
merged_table <- as.data.frame(merged_table)

# Normalize the dicermutationtype for each cancer's total mutation frequency
merged_table[1:5] <- lapply(merged_table[1:5], "*", merged_table$ratio)


merged_table$CancerType <- row.names(merged_table)
row.names(merged_table) <- NULL

### Melt the table by ratio/yvalue and cancer_type_acronynm

melted_table <- melt(
    data = merged_table, id.vars = c("ratio", "CancerType"),
    variable.name = "MutationType", value.name = "Height"
)

### Add a column to the melted table with number of samples (recycles the cancer
# type)

melted_table$CancerType_numSample <- paste(names(sample_count50),
    "(n=", sample_count50, ")",
    sep = ""
)


### Generate an index by decreasing ratio from the melted table
index <- order(melted_table$ratio, decreasing = TRUE)


#### Plot the data using index to sort x-axis(CANCERTYPE) based on sorted y-axis
# (% of total mutations)

plot_i <-
    qplot(
        x = factor(melted_table$CancerType_numSample[index],
            levels = unique(melted_table$CancerType_numSample[index])
        ),
        y = melted_table$Height[index],
        ylab = "",
        xlab = "",
        geom = "blank",
        fill = melted_table$MutationType[order(melted_table$ratio,
            decreasing = TRUE
        )]
    ) +
    geom_bar(stat = "identity") +
    theme_classic() +
    ylab("% Dicer Mutation(per tumor subtype)") +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(legend.position = c(0.70, 0.60), legend.title = element_blank()) +
    scale_y_continuous(breaks = c(0, 1, 2.5, 5, 10), limits = c(0, 12))


ggsave(
    filename = "team3/results/Figure_2B.pdf", plot = plot_i,
    width = 12, height = 8
)
