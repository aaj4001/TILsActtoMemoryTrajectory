library(monocle)
library(ggpubr)
library(ggplot2)
library(cowplot)

##################################################################################################################################
## Plotting Trajectory before Regression

load("../Monocle_Subpopulation.CD8CD3_Reproduced.rda")

Monocle_Subpopulation.CD8CD3@phenoData@data[["Response"]] = factor(Monocle_Subpopulation.CD8CD3@phenoData@data[["Response"]])
Monocle_Subpopulation.CD8CD3@phenoData@data[["Treatment_Status"]] = factor(Monocle_Subpopulation.CD8CD3@phenoData@data[["Treatment_Status"]],levels = c("Pre","Post"))

plot_cell_trajectory(Monocle_Subpopulation.CD8CD3, color_by = "State",cell_size = 0.75) + 
  facet_wrap(~Response+Treatment_Status, nrow = 1)+guides(colour = guide_legend(override.aes = list(size=5)))

Trajectories = plot_cell_trajectory(Monocle_Subpopulation.CD8CD3, color_by = "State",cell_size = 0.7) + 
  facet_wrap(~Response+Treatment_Status, nrow = 1)+guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(axis.title = element_blank(),axis.text = element_blank(),legend.position = "none",title = element_blank())

##################################################################################################################################
## Pie Chart before Regression

BeforeRegression = Monocle_Subpopulation.CD8CD3@phenoData@data[,c(1,2,3,13,14)]; colnames(BeforeRegression) = c("PatientID","Treatment_Status","Response","Before_Pseudo","Before_State")

PieCharts = data.frame(State = 1:5,
                       Before_NR_Pre = as.numeric(table(BeforeRegression$Before_State[BeforeRegression$Treatment_Status == "Pre" & BeforeRegression$Response == "Non-responder"])),
                       Before_NR_Post = as.numeric(table(BeforeRegression$Before_State[BeforeRegression$Treatment_Status == "Post" & BeforeRegression$Response == "Non-responder"])),
                       Before_R_Pre = as.numeric(table(BeforeRegression$Before_State[BeforeRegression$Treatment_Status == "Pre" & BeforeRegression$Response == "Responder"])),
                       Before_R_Post = as.numeric(table(BeforeRegression$Before_State[BeforeRegression$Treatment_Status == "Post" & BeforeRegression$Response == "Responder"])))
PieChartPercs = PieCharts
for(i in 2:5) PieChartPercs[,i] = round(PieChartPercs[,i]/sum(PieChartPercs[,i])*100)
print(PieChartPercs)

PieChart_Plots = list()
for(i in 2:5){
  PieChartRaw = data.frame(State = as.character(1:5),
                           Count = PieCharts[,i],
                           Perc = paste(round(PieCharts[,i]/sum(PieCharts[,i])*100),"%",sep = ""))
  # PieChart_Plots[[colnames(PieCharts)[i]]] = ggpie(PieChartRaw,"Count",label = "Perc",lab.pos = "in",color = "white",fill = "State",palette = get_palette("default",5)) +
  PieChart_Plots[[colnames(PieCharts)[i]]] = ggpie(PieChartRaw,"Count",label = "Perc",lab.pos = "in",color = "white",fill = "State",lab.font = c(0, "bold", "black"),palette = get_palette("default",5)) +
    theme(legend.position = "none") #+
  # ggtitle(colnames(PieCharts)[i])
}

plot_grid(plotlist = PieChart_Plots,nrow = 1)

Scale = 1.5
pdf("3AB_PlottingBeforeTrajectory_WithPie.pdf",height = 3.2*Scale, width = 5*Scale)
cowplot::plot_grid(Trajectories,
                   plot_grid(plotlist = PieChart_Plots,nrow = 1),nrow = 2,rel_heights = c(1.1,1))
dev.off()
