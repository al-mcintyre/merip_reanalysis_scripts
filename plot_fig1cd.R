library(ggplot2)
library(ComplexHeatmap)
library(circlize) #for colorRamp2

cluster <- TRUE
setwd('/Volumes/Markarian_501/merip_reanalysis/fig1b/')
peaks.overlap <- read.csv('peak_overlaps_summary.txt',sep=' ')
overlap.annotations <- read.csv('fig1cd_exp_summary.txt',sep=' ')
peaks.overlap$cellline1 <- sub("_", " ", peaks.overlap$cellline1)
peaks.overlap$cellline2 <- sub("_", " ", peaks.overlap$cellline2)

peaks.overlap$exp1 <- as.character(peaks.overlap$exp1)
peaks.overlap$exp2 <- as.character(peaks.overlap$exp2)
for (exp in c(1,2)){
  peaks.overlap[which(peaks.overlap[paste0("exp",exp)] == "Lichinchi" & 
                        peaks.overlap[paste0("cellline",exp)] == "MT CD4"),paste0("exp",exp)] <- "Lichinchi(a)"
  peaks.overlap[which(peaks.overlap[paste0("exp",exp)] == "Lichinchi" & 
                        peaks.overlap[paste0("cellline",exp)] == "HEK293T"),paste0("exp",exp)] <- "Lichinchi(b)"  
}

for (cellline in unique(peaks.overlap$cellline1)){
  cell.peak.ol <- peaks.overlap[which(peaks.overlap$cellline1 == cellline & peaks.overlap$cellline2 == cellline),]
  perc.overlap.mat <- matrix(data = NA,nrow=length(unique(cell.peak.ol$exp1)),ncol=length(unique(cell.peak.ol$exp2)))
  perc.overlap.mat <- as.data.frame(perc.overlap.mat)
  rownames(perc.overlap.mat) <- unique(cell.peak.ol$exp1)
  colnames(perc.overlap.mat) <- unique(cell.peak.ol$exp1)
  for (exp1 in cell.peak.ol$exp1){
    perc.overlap.mat[exp1,exp1] <- 100
    for (exp2 in cell.peak.ol[which(cell.peak.ol$exp1 == exp1),'exp2']){
      peak.overlap <- cell.peak.ol[which(cell.peak.ol$exp1 == exp1 & cell.peak.ol$exp2 == exp2),'peak_overlap']
      total.peaks <- cell.peak.ol[which(cell.peak.ol$exp1 == exp1 & cell.peak.ol$exp2 == exp2),'total_peaks']
      perc.overlap.mat[exp1,exp2] <- peak.overlap*100/total.peaks
    }
  }

  perc.overlap.sum <- as.vector(as.matrix(perc.overlap.mat))
  perc.overlap.sum <- perc.overlap.sum[which(perc.overlap.sum != 100)]
  summary <- paste(summary(perc.overlap.sum))
  writeLines(paste(cellline,'1st quartile = ',summary[2],', 3rd quartile = ',summary[5]))
  
    #could add the total #peaks per comparison
  # pdf(paste0("fig1c_",cellline,"_peak_heatmap.pdf"),width=3.5,height=3)
  # ch = Heatmap(perc.overlap.mat,col=colorRamp2(c(0,100),c("white","#89043d")), cluster_rows = FALSE, cluster_columns = FALSE, 
  #         column_names_side = "top",row_names_side = "left",name="% peaks\noverlapped",rect_gp = gpar(col = "white", lwd = 2),
  #         column_title = "Experiment 2",row_title = "Experiment 1",column_title_gp = gpar(fontsize=12),row_title_gp = gpar(fontsize=12)) 
  # draw(ch, column_title=cellline,column_title_gp = gpar(fontsize=14))
  # dev.off()
}


# for (species in unique(peaks.overlap$cellline1)){
#   cell.peak.ol <- peaks.overlap[which(peaks.overlap$species == species & peaks.overlap$species2 == cellline),]
#   perc.overlap.mat <- matrix(data = NA,nrow=length(unique(cell.peak.ol$exp1)),ncol=length(unique(cell.peak.ol$exp2)))
#   perc.overlap.mat <- as.data.frame(perc.overlap.mat)
#   rownames(perc.overlap.mat) <- unique(cell.peak.ol$exp1)
#   colnames(perc.overlap.mat) <- unique(cell.peak.ol$exp1)
#   for (exp1 in cell.peak.ol$exp1){
#     perc.overlap.mat[exp1,exp1] <- 100
#     for (exp2 in cell.peak.ol[which(cell.peak.ol$exp1 == exp1),'exp2']){
#       peak.overlap <- cell.peak.ol[which(cell.peak.ol$exp1 == exp1 & cell.peak.ol$exp2 == exp2),'peak_overlap']
#       total.peaks <- cell.peak.ol[which(cell.peak.ol$exp1 == exp1 & cell.peak.ol$exp2 == exp2),'total_peaks']
#       perc.overlap.mat[exp1,exp2] <- peak.overlap*100/total.peaks
#     }
#   }
#   
#   #could add the total #peaks per comparison
#   pdf(paste0("fig1d_",species,"_peak_heatmap.pdf"),width=3.5,height=3)
#   ch = Heatmap(perc.overlap.mat,col=colorRamp2(c(0,100),c("white","#89043d")), cluster_rows = FALSE, cluster_columns = FALSE, 
#                column_names_side = "top",row_names_side = "left",name="% peaks\noverlapped",rect_gp = gpar(col = "white", lwd = 2),
#                column_title = "Experiment 2",row_title = "Experiment 1",column_title_gp = gpar(fontsize=12),row_title_gp = gpar(fontsize=12)) 
#   draw(ch, column_title=cellline,column_title_gp = gpar(fontsize=14))
#   dev.off()
# }


for (species in c("hg38")){
  cell.lines <- sub("_"," ",overlap.annotations[which(overlap.annotations$species == species),"cell_line"])
  cell.peak.ol <- peaks.overlap[which(peaks.overlap$cellline1 %in% cell.lines),]
  cell.peak.ol$explabel1 <- paste0(cell.peak.ol$exp1," (",cell.peak.ol$cellline1,")")
  cell.peak.ol$explabel2 <- paste0(cell.peak.ol$exp2," (",cell.peak.ol$cellline2,")")
  perc.overlap.mat <- matrix(data = NA,nrow=length(unique(cell.peak.ol$explabel1)),ncol=length(unique(cell.peak.ol$explabel2)))
  perc.overlap.mat <- as.data.frame(perc.overlap.mat)
  rownames(perc.overlap.mat) <- unique(cell.peak.ol$explabel1) #paste(cell.peak.ol$exp1,cell.peak.ol$cellline1))
  colnames(perc.overlap.mat) <- unique(cell.peak.ol$explabel2) #paste(cell.peak.ol$exp2,cell.peak.ol$cellline2))
  for (exp1 in cell.peak.ol$explabel1){
    perc.overlap.mat[exp1,exp1] <- 100
    for (exp2 in cell.peak.ol[which(cell.peak.ol$explabel1 == exp1),'explabel2']){
      peak.overlap <- cell.peak.ol[which(cell.peak.ol$explabel1 == exp1 & cell.peak.ol$explabel2 == exp2),'peak_overlap']
      total.peaks <- cell.peak.ol[which(cell.peak.ol$explabel1 == exp1 & cell.peak.ol$explabel2 == exp2),'total_peaks']
      perc.overlap.mat[exp1,exp2] <- peak.overlap*100/total.peaks
    }
  }
}

overlap.annotations$explabel <- sub("_"," ",paste0(overlap.annotations$experiment_label," (",overlap.annotations$cell_line,")"))
row.order <- sort(overlap.annotations[which(overlap.annotations$explabel %in% rownames(perc.overlap.mat)),]$explabel)
perc.overlap.mat <- perc.overlap.mat[row.order,row.order]
o.a <- overlap.annotations[which(overlap.annotations$explabel %in% rownames(perc.overlap.mat)),]
row.names(o.a) <- o.a$explabel
o.a <- o.a[row.order,]
cl <- sort(unique(sub("_"," ",o.a$cell_line)))
cell.colours <- list("cell.line"= setNames(ggpubr::get_palette(c("#3f1a1c","#d78521","#e6af2e","#ede5a6"),length(cl)),cl))
st <- unique(o.a$experiment_label)
study.colours <- list("study"= setNames(ggpubr::get_palette("lancet",length(st)),st)) #c("#003f20","#bce5d2")

celllines = rowAnnotation(df = data.frame("cell.line" = sub("_"," ",o.a$cell_line)),
                              col = cell.colours) 
studies = rowAnnotation(df = data.frame("study" = o.a$experiment_label),col = study.colours) 

pdf(paste0("fig1d_clustered_peak_heatmap.pdf"),width=8,height=6)
ch = Heatmap(perc.overlap.mat,col=colorRamp2(c(0,100),c("white","#89043d")), cluster_rows = TRUE, cluster_columns = TRUE, 
             column_names_side = "top",row_names_side = "left",name="% peaks\noverlapped",rect_gp = gpar(col = "white", lwd = 2),
             column_title = "Experiment 2",row_title = "Experiment 1",column_title_gp = gpar(fontsize=12),row_title_gp = gpar(fontsize=12)) 
draw(ch + celllines + studies, heatmap_legend_side = "right")
dev.off()

perc.overlap.sum <- as.vector(as.matrix(perc.overlap.mat))
perc.overlap.sum <- perc.overlap.sum[which(perc.overlap.sum != 100)]
summary(perc.overlap.sum)

