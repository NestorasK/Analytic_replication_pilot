########## -----------GSEA analysis using clusterProfiler
########## -----------http://guangchuangyu.github.io/2015/05/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/
library("clusterProfiler")
require("DOSE")

########## -----------GSEA analysis of AR-targets in BML1 KD expression data.
# -pathway loading
gda <- read.delim("team4/data/AR_targets.txt")
AR2gene <- gda[, c("termName", "geneSymbol")]

############ -geneList loading
BMILfc <- read.delim("team4/data/BMILFC.rnk", header = FALSE)
head(BMILfc)
# Convert it to named Vector ##
BMILfc.vec <- BMILfc$V2
names(BMILfc.vec) <- BMILfc$V1
head(BMILfc.vec)

############ -running score and preranked list of GSEA result
fcBMIL <- GSEA(sort(BMILfc.vec, decreasing = TRUE),
    TERM2GENE = AR2gene
)
gseaplot(fcBMIL, "AR_Target")

############ -----BML1t.rnk
BML1t <- read.delim("team4/data/BML1t.rnk", header = FALSE)
head(BML1t)
# Convert it to named Vector ##
BML1t.vec <- BML1t$V2
names(BML1t.vec) <- BML1t$V1
head(BML1t.vec)

tBML1 = GSEA(sort(BML1t.vec, decreasing = TRUE),
    TERM2GENE = AR2gene
)
gseaplot(tBML1, "AR_Target")

############ ------------BML1p.rnk
BML1p <- read.delim("team4/data/BML1p.rnk", header = FALSE)
head(BML1p)
# Convert it to named Vector ##
BML1p.vec <- BML1p$V2
names(BML1p.vec) <- BML1p$V1
head(BML1p.vec)

pBML1 <- GSEA(sort(BML1p.vec, decreasing = TRUE), TERM2GENE = AR2gene)
gseaplot(pBML1, "AR_Target")

############ ------------BML1ap.rnk
BML1ap <- read.delim("team4/data/BML1ap.rnk", header = FALSE)
head(BML1ap)
# Convert it to named Vector ##
BML1ap.vec <- BML1ap$V2
names(BML1ap.vec) <- BML1ap$V1
head(BML1ap.vec)

apBML1 <- GSEA(sort(BML1ap.vec, decreasing = TRUE),
    TERM2GENE = AR2gene
)
gseaplot(apBML1, "AR_Target")


########### ----------GSEA analysis of AR-targets in RNF2fc KD expression data.
RNF2fc <- read.delim("team4/data/RNF2fc.rnk", header = FALSE)
head(RNF2fc)
RNF2fc.vec <- RNF2fc$V2
names(RNF2fc.vec) <- RNF2fc$V1
head(RNF2fc.vec)

fcRNF2 <- GSEA(sort(RNF2fc.vec, decreasing = TRUE),
    TERM2GENE = AR2gene
)
gseaplot(fcRNF2, "AR_Target")
