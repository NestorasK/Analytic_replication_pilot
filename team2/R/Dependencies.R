# Install EdgeR, statmod; RUN IF NEEDED
if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
  }
BiocManager::install("edgeR")
install.packages("statmod")
install.packages("devtools")
library(devtools, help, pos = 2, lib.loc = NULL)
install_github("vqv/ggbiplot")
install.packages("pvclust")
# install gplots for heatmap plot
install.packages("gplots")
install.packages("RColorBrewer")
install.packages("bcv")
