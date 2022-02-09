####################################################################################################
######### R script for bootstrap resampling analysis of TCGA and MSK-IMPACT cancer data ############
####################################################################################################

## This code is to perform a bootstrap resampling analysis for TCGA and MSK-IMPACT data for 
## enrichment of hotspot mutations. The script randomly resamples the number of cancer cases
## for every cancer type 10000 times to generate a distribution of observed hotspot mutations
## from the random draw. The random draw is performed on the total TCGA cancer cases without
## replacement.


rm(list=ls())

####### Load requried libraries #####

library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)


##### Load file #####

##### tcga is the file for TCGA cancer data ####

tcga <- read.table("TCGA_clinical_final.out", header=TRUE, quote=NULL)


#load table first then set the pvaltable
pvaltable <- rep(42, times = length(unique(tcga$CANCER_TYPE_ACRONYM)))



###Find POLE cases and replace the tcga/msk table with POLE cASES excluded for 2nd iteration of resampling

#pole<-c("COADREAD_POLE","UCEC_POLE","ESCA_POLE","STAD_POLE","ESCA_MSI","COADREAD_MSI","UCEC_MSI","STAD_MSI","")

#unique(tcga$CANCER_TYPE_ACRONYM)

#tcga<-tcga[-(tcga$CANCER_TYPE_ACRONYM %in% pole),]

##### msk is the file for MSK-IMPACT cancer data #### comment out the file as necessary,
##### can run only data from one file at a time

# msk <- read.table("MSK_clinical_final.out", header=TRUE, quote=NULL)

### Replace TCGA and MSK vice-versa in the following script to generate outputs for these data ###

perm = 10000

for (i in unique(tcga$CANCER_TYPE_ACRONYM)) {
  
  tempx=paste(i,sep="")
  
  type<- tcga[grep(tempx, tcga$CANCER_TYPE_ACRONYM), ]
  
  
  
  tcount <- nrow(type)
  multi <- type[grep("multi_with_hotspot", type$dicer_muts), ]
  mcount <- nrow(multi)
  hotspot <- type[grep("only_hotspot", type$dicer_muts), ]
  hcount <- nrow(hotspot)
  bial <- type[grep("bialelic", type$dicer_muts), ]
  bcount <- nrow(bial)
  sum <- mcount+hcount+bcount
  
  
  val <- capture.output(for (rand in 1:10000) {
    
    rand <- sample_n(tcga, tcount)
    mul <- rand[grep("multi_with_hotspot", rand$dicer_muts), ]
    mulcount <- nrow(mul)
    hot <- rand[grep("only_hotspot", rand$dicer_muts), ]
    hotcount <- nrow(hot)
    bi <- rand[grep("bialelic", rand$dicer_muts), ]
    bialcount <- nrow(bi)
    
    allcount <- mulcount+hotcount+bialcount
    cat(allcount,"\n", sep="")
    
  })
  
  val <- as.numeric(as.character(val))
  val <- melt(val)
  
  sum <- as.numeric(as.character(sum))
  sum <- melt(sum)
  
  sub <- subset(val, value > sum$value)
  
  rand <- nrow(sub)
  rand <- melt(rand)
  
  pval <- rand$value / perm
  
  top5 <- quantile(val$value, 0.95)
  top5 <- melt(top5)
  
    
 pval <- p.adjust(pval, method = "bonferroni", n=73)
  
 # pval <- formatC(pval, format = "e", digits = 2)
  
 pval <- formatC(pval, format = NULL, digits = 2)
  
  bot5 <- quantile(val$value, 0.05)
  bot5 <- melt(bot5)
  
  ggplot(val, aes(val$value))+
    geom_histogram(binwidth=.5, position="dodge") + theme_bw()+
    geom_vline(data=val, aes(xintercept = mean(val$value, na.rm=TRUE)),col='darkorange',size=2.5)+
    geom_vline(data=sum, aes(xintercept = sum$value),col='blue1',size=2.5)+
    geom_vline(data=sum, aes(xintercept = top5$value),col='red',size=2.5)+
    geom_vline(data=sum, aes(xintercept = bot5$value),col='black',size=2.5)+
    annotate(geom="text", x=0.5, y=280, label = paste("P=", sprintf(pval)))+
    xlab("Hotspot mutations")+
    ggtitle(paste(i, "(",tcount,")")) +
    theme(plot.title = element_text(size = 14))+
    theme(axis.text.x = element_text(size=16))+
    theme(axis.text.y = element_text(size=16))
  ggsave(filename=paste("TCGA_",i,".pdf",sep=""))
  
}
write.csv(x = pvaltable, file = "tcga_pval.csv")


##### END #####
