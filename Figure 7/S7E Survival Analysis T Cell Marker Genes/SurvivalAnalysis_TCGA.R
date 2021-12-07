## ---------------------------
##
## Script name: TCellStates_SurvivalAnalysis_TCGA_v3_KM.R
##
## Purpose of script: Survival Analysis of TCGA data - Ruth's method integrated
##                    Added function to query survival of individual marker genes
##
## Author: Abhi Jaiswal
##
## Date Created: 2021-09-01
##
## Email: ajaiswal1995@gmail.com
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------
##
## Packages:

library(readr)
library(readxl)
library(GSVA)
library(pheatmap)
library(ggpubr)
library(fgsea)
library(xCell)
library(survminer)
library(survival)
library(cowplot)
library(msigdb)

####################################################################################################################################
## Imports Annotation File


ClinicalInfo <- read_delim("../../Skin Cancer Clinical Information.txt",
                           "\t", escape_double = FALSE, trim_ws = TRUE)

ClinicalInfo_Survival = subset(ClinicalInfo,select = c("case_submitter_id",
                                                       "days_to_last_follow_up",
                                                       "days_to_death",
                                                       "vital_status",
                                                       "ajcc_pathologic_stage",
                                                       "treatment_type",
                                                       "age_at_index"))

colnames(ClinicalInfo_Survival) = c("ID","Time","TOD","LD","Stage","Tx","Age")

ClinicalInfo_Survival = ClinicalInfo_Survival[seq(1,nrow(ClinicalInfo_Survival),by = 2),]

Tumor_Pheno = read_delim("../../SkinCancer_gdc_sample_sheet.2021-04-14.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)[,c(6,8,5)]
Tumor_Pheno = Tumor_Pheno[order(Tumor_Pheno$`Case ID`),]
Tumor_Pheno = Tumor_Pheno[seq(1,nrow(Tumor_Pheno),by = 2),]
colnames(Tumor_Pheno) = c("CaseID","SampleType","ProjectID")

ClinicalInfo_Survival = merge(ClinicalInfo_Survival,Tumor_Pheno,by = 1,all.x = TRUE)

OfficialTime = numeric(nrow(ClinicalInfo_Survival))
i = 1
for(i in 1:length(OfficialTime)){
  if(ClinicalInfo_Survival$TOD[i]!="'--") {OfficialTime[i] = as.numeric(as.character(ClinicalInfo_Survival$TOD[i]))
  } else if (ClinicalInfo_Survival$TOD[i]=="'--"&ClinicalInfo_Survival$Time[i]!="'--") {OfficialTime[i] = as.numeric(as.character(ClinicalInfo_Survival$Time[i]))
  } else if (ClinicalInfo_Survival$TOD[i]=="'--"&ClinicalInfo_Survival$Time[i]=="'--") {OfficialTime[i] = NA}
}

ClinicalInfo_Survival$TimeFinal = OfficialTime

ClinicalInfo_Survival = ClinicalInfo_Survival[!is.na(ClinicalInfo_Survival$TimeFinal),c(-2,-3)]
ClinicalInfo_Survival$LD = ClinicalInfo_Survival$LD=="Dead"
ClinicalInfo_Survival$TimeFinal = ClinicalInfo_Survival$TimeFinal/30.417

ClinicalInfo_Survival = subset(ClinicalInfo_Survival,SampleType=="Metastatic")

####################################################################################################################################
## Imports FPKM File

FPKM <- read_csv("../../SkinCancerTCGA_FPKMAnnotated.csv")

FPKM_GSVA = data.matrix(FPKM[,-c(1,2)])

FPKM_GSVA = FPKM_GSVA[,colnames(FPKM_GSVA) %in% ClinicalInfo_Survival$ID]
rownames(FPKM_GSVA) = FPKM$hgnc_symbol

####################################################################################################################################
## GSVA of FPKM Files

# HacohenSignature = msigdb::read.gmt("../../TCellSigs_NE.gmt")$genesets[c(1,5,9,7,11,3)]
# YF = msigdb::read.gmt("YF_RNASeq_v2.gmt")$genesets[c(1,3)]
# 
# SignatureScored = YF

Pal = wesanderson::wes_palette(name = "Darjeeling1",3)[c(2,3,1)]


SortingQuantiles <- function(InputGSVAFile){
  QuantileSorted = data.frame(Sample = colnames(InputGSVAFile))
  for(i in rownames(InputGSVAFile)){
    
    GSVAScore = InputGSVAFile[i,]
    GSVAQuantile = rep("Low",length(GSVAScore));names(GSVAQuantile) = names(GSVAScore)
    GSVAQuantile[GSVAScore >= quantile(GSVAScore)[2]] = "Medium"
    GSVAQuantile[GSVAScore >= quantile(GSVAScore)[4]] = "High"
    
    GSVAQuantile = factor(GSVAQuantile,levels = c("Low","Medium","High"),ordered = TRUE)
    
    QuantileSorted[i] = GSVAQuantile
  }
  QuantileSorted
}

PlotSurvival <- function(SignatureScored,Clinical = ClinicalInfo_Survival,Vertical = FALSE){
  
  GSVAResult = gsva(FPKM_GSVA,SignatureScored,method = "ssgsea")
  Clinical = Clinical[Clinical$ID%in%colnames(GSVAResult),]
  
  GSVA_Quantile = SortingQuantiles(GSVAResult)
  
  Clinical_Quantile = merge(Clinical,GSVA_Quantile,by = 1)
  
  IndPlot = list()
  PVals = numeric()
  GroupSize = list()
  j = rownames(GSVAResult)[1]

  for(j in rownames(GSVAResult)){
    
    HighMedLow = Clinical_Quantile[[j]]
    HML = numeric()
    for(h in c("Low","Medium","High")){
      HML[h] = sum(HighMedLow==h)
    }
    GroupSize[[j]] = HML

    Formula = as.formula(paste0("Surv(time = TimeFinal,event = LD)~",j))
    Fit = do.call(survfit,list(formula = Formula,data = Clinical_Quantile))
    Diff = do.call(survdiff,list(formula = Formula,data = Clinical_Quantile))
    
    PVals[j]  <- pchisq(Diff$chisq, length(levels(HighMedLow))-1, lower.tail=FALSE)

    IndPlot[[j]] = ggsurvplot(Fit,font.tickslab = c(4, "plain", "black"), censor.size = 3,
                              palette = Pal,pval = FALSE,title = j,xlab = "Time (Months)",
                              legend = "none",pval.method = FALSE)$plot
  }
  
  PValsBH = p.adjust(PVals,method = "BH")
  print(PValsBH)
  print(gtools::stars.pval(PValsBH))
  
  cat("\n\n\n")
  
  print(sapply(GroupSize,function(x) x))
  
  if(Vertical){
    IndPlot[[nrow(GSVAResult)]] = IndPlot[[nrow(GSVAResult)]] + theme(axis.title = element_blank())
    for(k in 1:(nrow(GSVAResult)-1)) IndPlot[[k]] = IndPlot[[k]] + theme(axis.title = element_blank(),axis.text = element_blank())
    for(k in 1:(nrow(GSVAResult))) IndPlot[[k]] = IndPlot[[k]] + theme(plot.title = element_text(size = 6))
    
    CombinedPlot = plot_grid(plotlist = IndPlot,ncol = 1,align = "hv")
  } else{
    IndPlot[[1]] = IndPlot[[1]] + theme(axis.title = element_blank())
    for(k in 2:nrow(GSVAResult)) IndPlot[[k]] = IndPlot[[k]] + theme(axis.title = element_blank(),axis.text = element_blank())

    CombinedPlot = plot_grid(plotlist = IndPlot,nrow = 1,align = "hv")
  }
  CombinedPlot
}


PlotSurvival_Markers <- function(Markers,Clinical = ClinicalInfo_Survival,Vertical = FALSE){
  if(length(Markers)==1) {Exprs = t(data.matrix(FPKM_GSVA[Markers,])); rownames(Exprs) = Markers
  } else Exprs =   data.matrix(FPKM_GSVA[Markers,])

  Clinical = Clinical[Clinical$ID%in%colnames(Exprs),]
  
  Exprs_Quantile = SortingQuantiles(Exprs)
  
  Clinical_Quantile = merge(Clinical,Exprs_Quantile,by = 1)
  
  IndPlot = list()
  PVals = numeric()
  GroupSize = list()

  for(j in Markers){
    
    HighMedLow = Clinical_Quantile[[j]]
    HML = numeric()
    for(h in c("Low","Medium","High")){
      HML[h] = sum(HighMedLow==h)
    }
    GroupSize[[j]] = HML
    
    Formula = as.formula(paste0("Surv(time = TimeFinal,event = LD)~",j))
    Fit = do.call(survfit,list(formula = Formula,data = Clinical_Quantile))
    Diff = do.call(survdiff,list(formula = Formula,data = Clinical_Quantile))
    
    PVals[j]  <- pchisq(Diff$chisq, length(levels(HighMedLow))-1, lower.tail=FALSE)
    
    IndPlot[[j]] = ggsurvplot(Fit,font.tickslab = c(4, "plain", "black"), censor.size = 3,
                              palette = Pal,pval = FALSE,title = j,xlab = "Time (Months)",
                              legend = "none",pval.method = FALSE)$plot
  }
  
  PValsBH = p.adjust(PVals,method = "BH")
  print(PValsBH)
  print(gtools::stars.pval(PValsBH))
  
  cat("\n\n\n")
  
  print(sapply(GroupSize,function(x) x))
  
  if(Vertical){
    IndPlot[[nrow(Exprs)]] = IndPlot[[nrow(Exprs)]] + theme(axis.title = element_blank())
    for(k in 1:(nrow(Exprs)-1)) IndPlot[[k]] = IndPlot[[k]] + theme(axis.title = element_blank(),axis.text = element_blank())
    for(k in 1:(nrow(Exprs))) IndPlot[[k]] = IndPlot[[k]] + theme(plot.title = element_text(size = 6))
    
    CombinedPlot = plot_grid(plotlist = IndPlot,ncol = 1,align = "hv")
  } else{
    IndPlot[[1]] = IndPlot[[1]] + theme(axis.title = element_blank())
    for(k in 2:nrow(Exprs)) IndPlot[[k]] = IndPlot[[k]] + theme(axis.title = element_blank(),axis.text = element_blank())
    
    CombinedPlot = plot_grid(plotlist = IndPlot,nrow = 1,align = "hv")
  }
  CombinedPlot
}
  

####################################################################################################################################
## Survival Analysis

## T Cell States
Scale = 2
## TCF7 PD1 Expression
pdf(file = "S7E_CD8_GZMB_PD1_TCF7_IFNG Survival TCGA.pdf",width = 8*Scale,height = 1.3*Scale)
PlotSurvival_Markers(c("CD8A","GZMB","PDCD1","TCF7","IFNG"))
dev.off()








