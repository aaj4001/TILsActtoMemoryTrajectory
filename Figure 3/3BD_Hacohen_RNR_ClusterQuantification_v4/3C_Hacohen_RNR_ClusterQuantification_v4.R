library(openxlsx)
library(nlme)
library(emmeans)
library(ggplot2)

Hacohen_Populations = read.xlsx("../../Hacohen_TCellPopulations_ENTPD1_Tox_rankCutoff.xlsx",sheet = "RawData")

Hacohen_Populations_Perc = unique(Hacohen_Populations$PatientID)

Hacohen_Populations_Perc = data.frame(Patient = rep(Hacohen_Populations_Perc,4),
                                      RNR = rep(c("Responder","Non-responder"),each = 64),
                                      PrePost = rep(c("Pre","Post","Pre","Post"),each = 32),
                                      Cluster1 = NA,Cluster3 = NA,Cluster5 = NA)

for(i in 1:128){
  Hacohen_PopSize = Hacohen_Populations[Hacohen_Populations$PatientID==Hacohen_Populations_Perc$Patient[i] & Hacohen_Populations$Response==Hacohen_Populations_Perc$RNR[i] & Hacohen_Populations$Treatment_Status == Hacohen_Populations_Perc$PrePost[i],]
      
  if(nrow(Hacohen_PopSize)!=0){
    Hacohen_Populations_Perc$Cluster1[i] = sum(Hacohen_PopSize$State==1)/nrow(Hacohen_PopSize) * 100
    Hacohen_Populations_Perc$Cluster3[i] = sum(Hacohen_PopSize$State==3)/nrow(Hacohen_PopSize) * 100
    Hacohen_Populations_Perc$Cluster5[i] = sum(Hacohen_PopSize$State==5)/nrow(Hacohen_PopSize) * 100

    # Hacohen_Populations_Perc$Cluster1[i] = sum(Hacohen_PopSize$State==1)
    # Hacohen_Populations_Perc$Cluster3[i] = sum(Hacohen_PopSize$State==3)
    # Hacohen_Populations_Perc$Cluster5[i] = sum(Hacohen_PopSize$State==5)
  }
}

Hacohen_Populations_Perc = Hacohen_Populations_Perc[rowSums(is.na(Hacohen_Populations_Perc[,4:6]))==0,]
Hacohen_Populations_Percv2 = reshape2::melt(Hacohen_Populations_Perc,variable.name = c("Cluster"),value.name = "Percentage")

Clusters = c("Cluster1","Cluster3","Cluster5")

for(i in Clusters){
  Hacohen_LME = lme(Percentage~RNR*PrePost,random = ~1|Patient,data = Hacohen_Populations_Percv2[Hacohen_Populations_Percv2$Cluster==i,])
  em = emmeans(Hacohen_LME,"pairwise"~PrePost|RNR,adjust = "none")
  em2 = emmeans(Hacohen_LME,"pairwise"~RNR|PrePost,adjust = "none")
  plot(em)
  
  print(i)
  print(summary(em2))
}

openxlsx::write.xlsx(Hacohen_Populations_Percv2,file = "3C_Hacohen_RNR_ClusterQuantification_v2.xlsx")
