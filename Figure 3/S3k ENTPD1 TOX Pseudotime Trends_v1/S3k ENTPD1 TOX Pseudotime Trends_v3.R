## ---------------------------
##
## Script name: S3k ENTPD1 TOX Pseudotime Trends_v1.R
##
## Purpose of script: Repeat analysis done in Figure 3b/c on ENTPD1 TOX Populations Separately
##
## Author: Abhi Jaiswal
##
## Date Created: 2021-08-17
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

library(openxlsx)
library(nlme)
library(emmeans)

####################################################################################################################################

Hacohen_Populations = read.xlsx("../../Hacohen_TCellPopulations_ENTPD1_Tox_rankCutoff.xlsx",sheet = "RawData")

ENTPD1TOX_Status = c(DN = "DN",ENTPD1 = "ENTPD1",TOX = "TOX",DP = "DP")

RvsNR_ENTPD1 = lapply(ENTPD1TOX_Status,function(x){
  Hacohen_Populations_Perc = unique(Hacohen_Populations$PatientID)
  Hacohen_Populations_Perc = data.frame(Patient = rep(Hacohen_Populations_Perc,4),
                                        RNR = rep(c("Responder","Non-responder"),each = 64),
                                        PrePost = rep(c("Pre","Post","Pre","Post"),each = 32),
                                        Cluster1 = NA,Cluster3 = NA,Cluster5 = NA)
  i = 4
  for(i in 1:128){
    Hacohen_PopSize = subset(Hacohen_Populations,
                             PatientID == Hacohen_Populations_Perc$Patient[i] & 
                               Response == Hacohen_Populations_Perc$RNR[i] & 
                               Treatment_Status == Hacohen_Populations_Perc$PrePost[i] &
                               Population == x)
    
    
    if(nrow(Hacohen_PopSize)!=0){
      Hacohen_Populations_Perc$Cluster1[i] = sum(Hacohen_PopSize$State==1)/nrow(Hacohen_PopSize) * 100
      Hacohen_Populations_Perc$Cluster3[i] = sum(Hacohen_PopSize$State==3)/nrow(Hacohen_PopSize) * 100
      Hacohen_Populations_Perc$Cluster5[i] = sum(Hacohen_PopSize$State==5)/nrow(Hacohen_PopSize) * 100
    }
  }
  
  Hacohen_Populations_Perc = Hacohen_Populations_Perc[rowSums(is.na(Hacohen_Populations_Perc[,4:6]))==0,]
  Hacohen_Populations_Percv2 = reshape2::melt(Hacohen_Populations_Perc,variable.name = c("Cluster"),value.name = "Percentage")
  
  Hacohen_LME = lme(Percentage~RNR*Cluster,random = ~1|Patient,data = Hacohen_Populations_Percv2)
  EMMeans = summary(emmeans(Hacohen_LME,"pairwise"~RNR*Cluster,adjust = "none"))$emmeans
  
  print(x)
  print(summary(emmeans(Hacohen_LME,"pairwise"~RNR|Cluster,adjust = "none"))$contrasts)
  
  
  cat("\n\n\n")
  
  Hacohen_Populations_Percv2$ENTPD1_TOX = x
  EMMeans$ENTPD1_TOX = x
  list(RawData = Hacohen_Populations_Percv2,EMMeans = EMMeans)
})

RvsNR_ENTPD1_Concat = lapply(RvsNR_ENTPD1,function(x) x$RawData)
RvsNR_ENTPD1_Concat = do.call(rbind,RvsNR_ENTPD1_Concat)

Hacohen_LME = lme(Percentage~RNR*Cluster,random = ~1|Patient /ENTPD1_TOX,data = RvsNR_ENTPD1_Concat)
print(summary(emmeans(Hacohen_LME,"pairwise"~RNR|Cluster,adjust = "none")))

RvsNR_ENTPD1_Concat_EMMeans = summary(emmeans(Hacohen_LME,"pairwise"~RNR|Cluster,adjust = "none"))$emmeans

RvsNR_ENTPD1 = list(RawData = dplyr::bind_rows(lapply(RvsNR_ENTPD1,function(x) x$RawData)),
                    EMMEans = dplyr::bind_rows(lapply(RvsNR_ENTPD1,function(x) x$EMMeans)))

RvsNR_ENTPD1_Concat = list(RawData = RvsNR_ENTPD1_Concat,
                           EMMeans = RvsNR_ENTPD1_Concat_EMMeans)


openxlsx::write.xlsx(RvsNR_ENTPD1,"S3k ENTPD1 TOX Pseudotime Trends_v1.xlsx")
openxlsx::write.xlsx(RvsNR_ENTPD1_Concat,"S3k ENTPD1 TOX Pseudotime Trends_v3.xlsx")
