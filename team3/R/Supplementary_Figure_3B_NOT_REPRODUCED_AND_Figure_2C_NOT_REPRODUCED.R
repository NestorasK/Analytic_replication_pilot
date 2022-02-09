rm(list = ls())
library(dplyr)
library(ggplot2)
library(reshape)
setwd("PATH_TO_Analytic_replication_pilot")

# Breast Cancer
msk <- read.table("team3/data/msk_clinical_final.out",
  header = TRUE, quote = NULL
)

set.seed(42) # set seed for the for loop

perm <- 10000 # number of permutations to execute


# cancer type acronyms were mislabeled in the tabble originally and so when
# added to this vector they didn't capture the correct types
cancer_types <- c(
  "breast_cancer", "endometrial_cancer",
  "small_cell_lung_cancer", "nsc_lung_cancer"
)
cancer_types <- as.factor(x = cancer_types)

ctypepval <- data.frame(
  Cancer_Types = rep(NA, times = length(cancer_types)),
  pval = rep(NA, times = length(cancer_types))
)

counter <- 0

for (i in unique(cancer_types)) {
  # selecting what is to be measured
  type <- msk[grepl(pattern = i, x = msk$CANCER_TYPE_ACRONYM), ]

  # setting the variables which need to be measured
  tcount <- nrow(type)
  multi <- type[grep("multi_with_hotspot", type$dicer_muts), ]
  mcount <- nrow(multi) # counting the types of multihotspot by dicer mutants
  hotspot <- type[grep("only_hotspot", type$dicer_muts), ]
  hcount <- nrow(hotspot)
  bial <- type[grep("bialelic", type$dicer_muts), ]
  bcount <- nrow(bial)
  sum <- mcount + hcount + bcount


  val <- capture.output(for (rand in 1:perm) {
    rand <- sample_n(msk, tcount)
    mul <- rand[grep("multi_with_hotspot", rand$dicer_muts), ]
    mulcount <- nrow(mul)
    hot <- rand[grep("only_hotspot", rand$dicer_muts), ]
    hotcount <- nrow(hot)
    bi <- rand[grep("bialelic", rand$dicer_muts), ]
    bialcount <- nrow(bi)

    allcount <- mulcount + hotcount + bialcount
    cat(allcount, "\n", sep = "")
  })

  # makes a numeric vector of the above for loop
  val <- as.numeric(as.character(val))

  val <- melt(val) # reformats the vector into a table

  sum <- as.numeric(as.character(sum))
  sum <- melt(sum)

  sub <- subset(val, value > sum$value)

  rand <- nrow(sub)
  rand <- melt(rand)


  pval <- rand$value / perm

  top5 <- quantile(val$value, 0.95)
  top5 <- melt(top5)


  pval <- p.adjust(pval, method = "bonferroni", n = 73)

  pval <- formatC(pval, format = NULL, digits = 2)

  counter <- counter + 1

  ctypepval$Cancer_Types[counter] <- i
  ctypepval$pval[counter] <- pval


  bot5 <- quantile(val$value, 0.05)
  bot5 <- melt(bot5)

  # plotting the data set after forming
  plot_i <-
    ggplot(val, aes(val$value)) +
    geom_histogram(binwidth = .5, position = "dodge") +
    theme_bw() +
    geom_vline(
      data = val,
      aes(xintercept = mean(val$value, na.rm = TRUE)),
      col = "darkorange", size = 2.5
    ) +
    geom_vline(
      data = sum, aes(xintercept = sum$value),
      col = "blue1", size = 2.5
    ) +
    geom_vline(
      data = sum, aes(xintercept = top5$value),
      col = "red", size = 2.5
    ) +
    geom_vline(
      data = sum, aes(xintercept = bot5$value),
      col = "black", size = 2.5
    ) +
    annotate(
      geom = "text", x = 0.5, y = 280,
      label = paste("P=", sprintf(pval))
    ) +
    xlab("Hotspot mutations") +
    ggtitle(label = paste(i, "(", tcount, ")")) +
    theme(plot.title = element_text(size = 14)) +
    theme(axis.text.x = element_text(size = 16)) +
    theme(axis.text.y = element_text(size = 16))

  ggsave(
    filename = paste(
      "team3/results/Supplementary_Figure_3B_",
      i, ".pdf",
      sep = ""
    ),
    width = 5, height = 5, plot = plot_i
  )
}
