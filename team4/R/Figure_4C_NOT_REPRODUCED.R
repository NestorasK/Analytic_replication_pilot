
## ----message=FALSE, warning=FALSE, include=FALSE-------------------------
rm(list = ls())
library("SummarizedExperiment")
library("dplyr")
library("TCGAbiolinks")
library("grid")

query.exp <- GDCquery(
  project = "TCGA-PRAD",
  legacy = TRUE,
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq",
  file.type = "results",
  experimental.strategy = "RNA-Seq"
)
GDCdownload(query.exp, directory = "team4/data/GDCdata")
prad.exp <- GDCprepare(
  query = query.exp, save = TRUE, directory = "team4/data/GDCdata",
  save.filename = "team4/data/GDCdata/pradExp.rda"
)

# get clinical data
dataClin <- GDCquery_clinic(project = "TCGA-PRAD", "clinical")
head(dataClin)
dim(dataClin)


## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
dataPrep <- TCGAanalyze_Preprocessing(object = prad.exp, cor.cut = 0.6)
dim(dataPrep)

dataNorm <- TCGAanalyze_Normalization(
  tabDF = dataPrep,
  geneInfo = geneInfo,
  method = "gcContent"
)
dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataNorm,
  method = "quantile",
  qnt.cut = 0.25
)
dim(dataFilt)

datFilt <- dataNorm %>%
  TCGAanalyze_Filtering(method = "varFilter") %>%
  TCGAanalyze_Filtering(method = "filter1") %>%
  TCGAanalyze_Filtering(method = "filter2", foldChange = 0.2)

data_Hc2 <- TCGAanalyze_Clustering(
  tabDF = datFilt,
  method = "consensus",
  methodHC = "ward.D2"
)

## ---------Differentially expression analysis (DEA) using edgeR package.

dataDEGs <- TCGAanalyze_DEA(
  mat1 = dataFilt[, dataFilt["BMI1", ] <
    mean(dataFilt["BMI1", ])],
  mat2 = dataFilt[, dataFilt["BMI1", ] >= mean(dataFilt["BMI1", ])],
  Cond1type = "Low",
  Cond2type = "High",
  fdr.cut = 0.01,
  logFC.cut = 1,
  method = "glmLRT"
)


############ ---------Add  cluster information to Summarized Experiment.
colData(prad.exp)$BMI1Exp <- rep(x = NA, times = ncol(dataFilt))
for (i in 1:ncol(dataFilt)) {
  if (dataFilt["BMI1", i] >= mean(dataFilt["BMI1", ])) {
    colData(prad.exp)$BMI1Exp[i] <- "High"
  } else {
    colData(prad.exp)$BMI1Exp[i] <- "Low"
  }
}


############ ---------Survival
TCGAanalyze_survival(
  data = colData(prad.exp),
  clusterCol = "BMI1Exp",
  main = "TCGA kaplan meier survival plot from consensus cluster",
  legend = "RNA Group", height = 10,
  risk.table = T, conf.int = F,
  color = c("black", "red", "blue", "green3"),
  filename = "team4/results/survival_BML1_expression_subtypes.png"
)

############ ---------
library("png")
library("grid")

TCGAvisualize_Heatmap(t(datFilt),
  col.metadata = colData(prad.exp)[, c("barcode", "BMI1Exp")],
  col.colors = list(BMI1Exp = c("High" = "black", "Low" = "red")),
  sortCol = "BMI1Exp",
  type = "expression", # sets default color
  scale = "row", # use z-scores for better visualization. Center gene expression level around 0.
  title = "Heatmap from concensus cluster",
  filename = "team4/results/BML1z_Heatmap.png",
  cluster_rows = TRUE,
  color.levels = colorRampPalette(c(
    "blue", "black",
    "yellow"
  ))(n = 11),
  extremes = seq(-5, 5, 1),
  cluster_columns = FALSE,
  width = 1000,
  height = 1000
)
