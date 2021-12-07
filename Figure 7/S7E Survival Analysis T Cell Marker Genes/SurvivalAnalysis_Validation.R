## ---------------------------
##
## Script name: SurvivalAnalysis_Validation_RuthIntegrated_v2.R
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
##    - AN ERROR WAS FOUND IN THIS CODE WHEREIN THE INPUT DATA MATRIX (TPM) was stored as a character matrix - this has been corrected, and TPM is now a numeric matrix
##   
##
## ---------------------------
##
## Packages:


####################################################################################################################################

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

################################################################################################################################################################
## Imports Datasets File

TPM <- data.frame(read_delim("../../aPD1_Melanoma_TPM.txt", "\t",escape_double = FALSE, trim_ws = TRUE))

## Filters on PD, CR, PR and assigns them to NR/R status
Annotation <- read_excel("../../aPD1_Melanoma_Annotation.xlsx")
Annotation_Filtered = merge(data.frame(TPM$X1),Annotation,by = 1)
Annotation_Filtered$BR = factor(Annotation_Filtered$BR, levels = c("PD", "PR", "MR", "CR", "SD"))
Annotation_Filtered$priorCTLA4 = as.logical(Annotation_Filtered$priorCTLA4)


Annotation_Filtered = subset(Annotation_Filtered,select = c("TPM.X1","dead","BR","PFS","OS","TimeToBR"))
Annotation_Filtered$dead = Annotation_Filtered$dead==1
colnames(Annotation_Filtered) = c("ID","LD","BR","PFS","OS","TimetoBR")
Annotation_Filtered$TimeFinal = Annotation_Filtered$OS/30.417

ClinicalInfo_Survival = Annotation_Filtered
################################################################################################################################################################
## Imports FPKM File

rownames(TPM) = TPM$X1; TPM = TPM[,-1];TPM = data.matrix(TPM)
FPKM_GSVA = t(TPM)
FPKM_GSVA = FPKM_GSVA[,colnames(FPKM_GSVA) %in% ClinicalInfo_Survival$ID]

################################################################################################################################################################
## GSVA of FPKM Files

# TCellSigs = msigdb::read.gmt("../../TCellSigs_NE.gmt")$genesets[c(1,5,9,7,11,3)]
# YF = msigdb::read.gmt("YF_RNASeq_v3.gmt")$genesets[c(1,3)]
# 
# SignatureScored = YF


Pal = wesanderson::wes_palette(name = "Darjeeling1",3)[c(2,3,1)]

SortingQuantiles <- function(InputGSVAFile,Percentage = 0.25){
  QuantileSorted = data.frame(Sample = colnames(InputGSVAFile))
  
  for(i in 1:nrow(InputGSVAFile)){
    Quantiles = InputGSVAFile[i,]
    Quantiles = sort(Quantiles)
    
    Low = 1:floor(Percentage*length(Quantiles))
    Low = data.frame(Patient = names(Quantiles)[Low],Status = "Low")
    
    High = ceiling((1-Percentage)*length(Quantiles)):length(Quantiles)
    High = data.frame(Patient = names(Quantiles)[High],Status = "High")
    
    Medium = (floor(Percentage*length(Quantiles))+1):(ceiling((1-Percentage)*length(Quantiles))-1)
    Medium = data.frame(Patient = names(Quantiles)[Medium],Status = "Medium")
    
    
    Annotated_sub = rbind(Low,Medium,High); colnames(Annotated_sub) = c("Sample",rownames(InputGSVAFile)[i])
    Annotated_sub[,2] = factor(Annotated_sub[,2],levels = c("Low","Medium","High"),ordered = TRUE)
    QuantileSorted = merge(QuantileSorted,Annotated_sub,by = 1)
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


################################################################################################################################################################
## Survival Analysis

## T Cell States
Scale = 2
## TCF7 PD1 Expression
pdf(file = "S7E_CD8_GZMB_PD1_TCF7_IFNG Survival Validation.pdf",width = 8*Scale,height = 1.3*Scale)
PlotSurvival_Markers(c("CD8A","GZMB","PDCD1","TCF7","IFNG"))
dev.off()




