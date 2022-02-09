rm(list = ls())
library(dplyr)
library(ggplot2)
library(reshape)
setwd("PATH_TO_Analytic_replication_pilot")

# Load file

msk <- read.table("team3/data/MSK_clinical_final.out",
  header = TRUE, quote = NULL
)


set.seed(42) # set seed for the for loop

perm <- 10000 # doing 10000 permiatations, could be higher or lower

MSK_pval <- data.frame(
  CANCER_TYPE_ACRONYM = rep(NA,
    times = length(unique(msk$CANCER_TYPE_ACRONYM))
  ),
  pval = rep(NA, times = length(unique(msk$CANCER_TYPE_ACRONYM)))
)


counter <- 0

for (i in unique(msk$CANCER_TYPE_ACRONYM)) {
  tempx <- paste(i, sep = "") # indexing the future variables

  # selecting what is to be measured
  type <- msk[grep(tempx, msk$CANCER_TYPE_ACRONYM), ]


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

  MSK_pval$CANCER_TYPE_ACRONYM[counter] <- i
  MSK_pval$pval[counter] <- pval

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
    ggtitle(paste(i, "(", tcount, ")")) +
    theme(plot.title = element_text(size = 14)) +
    theme(axis.text.x = element_text(size = 16)) +
    theme(axis.text.y = element_text(size = 16))

  # output is individual parts to supplemental figure 2
  ggsave(filename = paste("team3/results/Supplementary_Figure3a_", i,
    ".pdf",
    sep = ""
  ), width = 5, height = 5)
}
write.csv(x = MSK_pval, file = "team3/results/Supplementary_Figure3a_pval.csv")
