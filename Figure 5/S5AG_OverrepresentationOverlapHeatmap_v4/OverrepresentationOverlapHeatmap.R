library(Seurat)
library(msigdb)
library(pheatmap)
library(cowplot)

HacohenPseudo = readRDS("../../Monocle_CD8CD3.seurat.rds")

Hacohen_Exprs = data.matrix(GetAssayData(HacohenPseudo))
Hacohen_States = data.frame(row.names = rownames(HacohenPseudo@meta.data),
                            Response = HacohenPseudo@meta.data$Response,
                            Cluster = HacohenPseudo@meta.data$State)

################################################################################
## AJ Heatmap function
# Hacohen_DEG = Hacohen_DEG_C1
Hacohen_Heatmap <- function(Hacohen_DEG,RowstoPlot = NULL,Breaks = c(-2,2),
                            RowOrderManual = FALSE){
  
  Hacohen_Exprs_Genes = Hacohen_Exprs[Hacohen_DEG,]
  Hacohen_Exprs_Genes = Hacohen_Exprs_Genes[,order(Hacohen_States$Cluster,Hacohen_States$Response)]
  for(i in 1:nrow(Hacohen_Exprs_Genes)) Hacohen_Exprs_Genes[i,] = scale(Hacohen_Exprs_Genes[i,])
  Hacohen_Exprs_Genes[Hacohen_Exprs_Genes > Breaks[2]] = Breaks[2]
  Hacohen_Exprs_Genes[Hacohen_Exprs_Genes < Breaks[1]] = Breaks[1]
  
  HacohenStates_Plot = Hacohen_States[order(Hacohen_States$Cluster,Hacohen_States$Response),]
  
  # RowstoPlot = rownames(Hacohen_Exprs_Genes)[sample(1:nrow(Hacohen_Exprs_Genes),size = 20)]
  
  if(!is.null(RowstoPlot)) {rownames(Hacohen_Exprs_Genes)[!rownames(Hacohen_Exprs_Genes)%in%RowstoPlot] = ""
  } else{
    rownames(Hacohen_Exprs_Genes) = rep("",nrow(Hacohen_Exprs_Genes))
  }
  
  HMap_Colors = list(Cluster = scales::hue_pal()(5),
                     Response = c("turquoise","maroon2"))
  names(HMap_Colors$Cluster) = levels(HacohenStates_Plot$Cluster)
  names(HMap_Colors$Response) = levels(HacohenStates_Plot$Response)
  
  
  
  if(RowOrderManual){
    HMap = pheatmap(Hacohen_Exprs_Genes,scale = "none",cluster_cols = FALSE,show_colnames = FALSE,cluster_rows = FALSE,
                  silent = TRUE,annotation_col = HacohenStates_Plot,annotation_colors = HMap_Colors,
                  breaks = seq(Breaks[1],Breaks[2],length.out = 101))$gtable
  } else{
    HMap = pheatmap(Hacohen_Exprs_Genes,scale = "none",cluster_cols = FALSE,show_colnames = FALSE,
                    silent = TRUE,annotation_col = HacohenStates_Plot,annotation_colors = HMap_Colors,
                    breaks = seq(Breaks[1],Breaks[2],length.out = 101))$gtable
  }
  
  HMap
}

################################################################################
## Cluster 1 Overlap to effector

Cluster1_Up = read.gmt("../../Hacohen_StateDEGs_Human.gmt")$genesets[[1]]
QueryGenes = list(Tact_Exc = read.gmt("../../TCellSigs_Exc.gmt")$genesets[[3]],
                  YF_Eff = read.gmt("../../YellowFever_EffectorMemory_Exc.gmt")$genesets[[1]],
                  YF_EFF_RNA = read.gmt("../../YellowFever_EffectorMemory_RNASeq_v4Exc.gmt")$genesets[[1]])

QueryGenes = lapply(QueryGenes,function(x) intersect(Cluster1_Up,x))
QueryIntersects = list(Act_YF = intersect(QueryGenes[[1]],QueryGenes[[2]]),
                       Act_RNA = intersect(QueryGenes[[1]],QueryGenes[[3]]),
                       YF_RNA = intersect(QueryGenes[[2]],QueryGenes[[3]]))
  
HMap_Genes = list(Act = setdiff(QueryGenes[[1]],c(QueryGenes[[2]],QueryGenes[[3]])),
                  YF = setdiff(QueryGenes[[2]],c(QueryGenes[[1]],QueryGenes[[3]])),
                  RNA = setdiff(QueryGenes[[3]],c(QueryGenes[[2]],QueryGenes[[1]])),
                  Act_YF = setdiff(QueryIntersects[[1]],c(QueryIntersects[[2]],QueryIntersects[[3]])),
                  Act_RNA = setdiff(QueryIntersects[[2]],c(QueryIntersects[[1]],QueryIntersects[[3]])),
                  YF_RNA = setdiff(QueryIntersects[[3]],c(QueryIntersects[[1]],QueryIntersects[[2]])),
                  Int = intersect(intersect(QueryIntersects[[1]],QueryIntersects[[2]]),QueryIntersects[[3]]))
                  
HMap_Genes = do.call(c,HMap_Genes)

Scale = 2.5
png("OverrepresentationOverlapHeatmap_Cluster1Eff.png",width = 2.3*Scale, height = 4.2*Scale,res = 500, units = "in")
ggdraw(Hacohen_Heatmap(HMap_Genes,RowstoPlot = HMap_Genes,RowOrderManual = TRUE,Breaks = c(-1.5,1.5)))
dev.off()

################################################################################
## Cluster 3 Overlap to TRM

Cluster3_Up = read.gmt("../../Hacohen_StateDEGs_Human.gmt")$genesets[[3]]
TRM_MemoryExc = c(read.gmt("../../TCellSigs_Exc.gmt")$genesets[c(7,1)],
                  read.gmt("../../YellowFever_EffectorMemory_RNASeq_v4Exc.gmt")$genesets[3])
names(TRM_MemoryExc) = c("TRM","Mem","YF_RNA")

TRM_MemoryExc = lapply(TRM_MemoryExc,function(x) intersect(x,Cluster3_Up))
HMap_Genes = list(TRM = setdiff(TRM_MemoryExc[[1]],c(TRM_MemoryExc[[2]],TRM_MemoryExc[[3]])),
                  Mem = setdiff(TRM_MemoryExc[[2]],c(TRM_MemoryExc[[1]],TRM_MemoryExc[[3]])),
                  RNA = setdiff(TRM_MemoryExc[[3]],c(TRM_MemoryExc[[1]],TRM_MemoryExc[[2]])),
                  TRM_RNA = intersect(TRM_MemoryExc[[1]],TRM_MemoryExc[[3]]),
                  Mem_RNA = intersect(TRM_MemoryExc[[2]],TRM_MemoryExc[[3]]))
HMap_Genes = do.call(c,HMap_Genes)

Scale = 2.5
png("OverrepresentationOverlapHeatmap_Cluster3 TRM Mem.png",width = 2.3*Scale, height = 3*Scale,res = 500, units = "in")
ggdraw(Hacohen_Heatmap(HMap_Genes,RowstoPlot = HMap_Genes,RowOrderManual = TRUE,Breaks = c(-1.5,1.5)))
dev.off()

################################################################################
## Cluster Overlaps to Mem

ClusterMem = read.gmt("../../Hacohen_StateDEGs_Human.gmt")$genesets[1]
TMem = read.gmt("../../TCellSigs_Exc.gmt")$genesets[[1]]

ClusterMem = lapply(ClusterMem,function(x) intersect(x,TMem))
ClusterMem = unlist(ClusterMem)

Scale = 2.5
png("OverrepresentationOverlapHeatmap_Cluster1Memory.png",width = 3.6*Scale, height = 1*Scale,res = 500, units = "in")
ggdraw(Hacohen_Heatmap(ClusterMem,RowstoPlot = ClusterMem,RowOrderManual = TRUE))
dev.off()





