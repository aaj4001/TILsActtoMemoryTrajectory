library(Seurat)
library(msigdb)
library(pheatmap)
library(cowplot)
library(biomaRt)

HacohenGenes = read.gmt("../../Hacohen_StateDEGs_Human.gmt")$genesets
HacohenPseudo = readRDS("../../Monocle_CD8CD3.seurat.rds")

Hacohen_DEG_C1 = unique(unlist(HacohenGenes[1:2]))
Hacohen_DEG_C3 = unique(unlist(HacohenGenes[3:4]))
Hacohen_DEG_C5 = unique(unlist(HacohenGenes[5:6]))

Hacohen_Exprs = data.matrix(GetAssayData(HacohenPseudo))
Hacohen_States = data.frame(row.names = rownames(HacohenPseudo@meta.data),
                            Response = HacohenPseudo@meta.data$Response,
                            Cluster = HacohenPseudo@meta.data$State)

DEG_List = list(Hacohen_DEG_C1 = Hacohen_DEG_C1,
                Hacohen_DEG_C3 = Hacohen_DEG_C3,
                Hacohen_DEG_C5 = Hacohen_DEG_C5)

################################################################################
## AJ Heatmap function
# Hacohen_DEG_List = DEG_List;Breaks = c(-2,2);palette = NULL
Hacohen_Heatmap <- function(Hacohen_DEG_List,Breaks = c(-2,2),palette = NULL){
  DEG_Clustered = lapply(Hacohen_DEG_List,function(x){
    Hacohen_Exprs_Genes = Hacohen_Exprs[x,]
    Hacohen_Exprs_Genes = Hacohen_Exprs_Genes[,order(Hacohen_States$Cluster,Hacohen_States$Response)]
    for(i in 1:nrow(Hacohen_Exprs_Genes)) Hacohen_Exprs_Genes[i,] = scale(Hacohen_Exprs_Genes[i,])
    Hacohen_Exprs_Genes[Hacohen_Exprs_Genes > Breaks[2]] = Breaks[2]
    Hacohen_Exprs_Genes[Hacohen_Exprs_Genes < Breaks[1]] = Breaks[1]
    
    HacohenStates_Plot = Hacohen_States[order(Hacohen_States$Cluster,Hacohen_States$Response),]
    
    HMap = pheatmap(Hacohen_Exprs_Genes,scale = "none",cluster_cols = FALSE,show_colnames = FALSE,
                    silent = TRUE)
    x[HMap[["tree_row"]][["order"]]]
  })
  
  ## Method 1 - Make into 1 heatmap
  DEG_Clustered_v2 = data.frame(ClusterRow = unlist(sapply(1:length(Hacohen_DEG_List),function(x) rep(names(Hacohen_DEG_List)[x],length(Hacohen_DEG_List[[x]])))))
  rownames(DEG_Clustered_v2) = c(1:nrow(DEG_Clustered_v2))
  DEG_Clustered_v2$ClusterRow = factor(DEG_Clustered_v2$ClusterRow,levels = names(Hacohen_DEG_List))
  
  Exprs_Clustered = sapply(DEG_Clustered,function(x){
    Hacohen_Exprs_Genes = Hacohen_Exprs[x,]
    Hacohen_Exprs_Genes = Hacohen_Exprs_Genes[,order(Hacohen_States$Cluster,Hacohen_States$Response)]
    for(i in 1:nrow(Hacohen_Exprs_Genes)) Hacohen_Exprs_Genes[i,] = scale(Hacohen_Exprs_Genes[i,])
    Hacohen_Exprs_Genes[Hacohen_Exprs_Genes > Breaks[2]] = Breaks[2]
    Hacohen_Exprs_Genes[Hacohen_Exprs_Genes < Breaks[1]] = Breaks[1]
    
    Hacohen_Exprs_Genes
  })
  
  HMap_Combined = Exprs_Clustered[[1]]
  for(i in 2:length(Exprs_Clustered)) HMap_Combined = rbind(HMap_Combined,Exprs_Clustered[[i]])
  rownames(HMap_Combined) = 1:nrow(HMap_Combined)
  
  HacohenStates_Plot = Hacohen_States[order(Hacohen_States$Cluster,Hacohen_States$Response),]
  
  if(is.null(palette)) palette = rep("black",length(Hacohen_DEG_List))
  
  Gaps = sapply(DEG_Clustered,function(x) length(x))
  for(i in 2:length(Gaps)) Gaps[i] = Gaps[i] + Gaps[(i-1)]
  Gaps = Gaps[-length(Gaps)] 
  
  HMap_Colors = list(Cluster = scales::hue_pal()(5),
                     Response = c("turquoise","maroon2"),
                     ClusterRow = palette)
  
  names(HMap_Colors$Cluster) = levels(HacohenStates_Plot$Cluster)
  names(HMap_Colors$Response) = levels(HacohenStates_Plot$Response)
  names(HMap_Colors$ClusterRow) = levels(DEG_Clustered_v2$ClusterRow)
  
  HMap = pheatmap(HMap_Combined,scale = "none",cluster_cols = FALSE,cluster_rows = FALSE,show_colnames = FALSE,show_rownames = FALSE,gaps_row = Gaps,
                  silent = TRUE,annotation_col = HacohenStates_Plot,annotation_row = DEG_Clustered_v2,annotation_colors = HMap_Colors)$gtable
  HMap
}

pal = scales::hue_pal()(5)[c(1,3,5)]
HMap_Combined = ggdraw(Hacohen_Heatmap(DEG_List,palette = pal))                                                               

Scale = 2
png("Hacohen_HMap_v3.png",width = 3*Scale, height = 5*Scale,res = 500, units = "in")
HMap_Combined
dev.off()


########################################################################################################################
## Sorts top 25 DEGs up and down for each cluster
HacohenPseudo = readRDS("../../Monocle_CD8CD3.seurat.rds")

Cluster1DEG = FindMarkers(HacohenPseudo,ident.1 = "1", ident.2 = c("3","5"), group.by = "State",min.pct = 0.25)
Cluster3DEG = FindMarkers(HacohenPseudo,ident.1 = "3", ident.2 = c("1","5"), group.by = "State",min.pct = 0.25)
Cluster5DEG = FindMarkers(HacohenPseudo,ident.1 = "5", ident.2 = c("3","1"), group.by = "State",min.pct = 0.25)

FDRCutoff = 0.01
Cluster1 = list(Up = Cluster1DEG[Cluster1DEG$avg_logFC>0 & Cluster1DEG$p_val_adj<FDRCutoff,],
                Dn = Cluster1DEG[Cluster1DEG$avg_logFC<0 & Cluster1DEG$p_val_adj<FDRCutoff,])
Cluster3 = list(Up = Cluster3DEG[Cluster3DEG$avg_logFC>0 & Cluster3DEG$p_val_adj<FDRCutoff,],
                Dn = Cluster3DEG[Cluster3DEG$avg_logFC<0 & Cluster3DEG$p_val_adj<FDRCutoff,])
Cluster5 = list(Up = Cluster5DEG[Cluster5DEG$avg_logFC>0 & Cluster5DEG$p_val_adj<FDRCutoff,],
                Dn = Cluster5DEG[Cluster5DEG$avg_logFC<0 & Cluster5DEG$p_val_adj<FDRCutoff,])

## Biomart protein coding filter
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

Cluster1 = lapply(Cluster1,function(x){
  ProteinCoding = getBM(mart = mart,attributes = c("external_gene_name","gene_biotype"),
                 filters = "external_gene_name",values = rownames(x))
  merge(ProteinCoding,x,by.x = 1,by.y = 0)
})
Cluster3 = lapply(Cluster3,function(x){
  ProteinCoding = getBM(mart = mart,attributes = c("external_gene_name","gene_biotype"),
                        filters = "external_gene_name",values = rownames(x))
  merge(ProteinCoding,x,by.x = 1,by.y = 0)
})
Cluster5 = lapply(Cluster5,function(x){
  ProteinCoding = getBM(mart = mart,attributes = c("external_gene_name","gene_biotype"),
                        filters = "external_gene_name",values = rownames(x))
  merge(ProteinCoding,x,by.x = 1,by.y = 0)
})

DEGClusters = list(Cluster1 = Cluster1,Cluster3 = Cluster3,Cluster5 = Cluster5)
DEGClusters = unlist(DEGClusters,recursive = FALSE)

DEGClusters = sapply(DEGClusters,function(x){
  x = subset(x,gene_biotype == "protein_coding")
  List = x$external_gene_name[order(abs(x$avg_logFC),decreasing = TRUE)]
  List[1:20]
})
