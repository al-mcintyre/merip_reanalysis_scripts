library(ggplot2)
library(ComplexHeatmap)
library(circlize) #for colorRamp2
#also requires ggpubr

args <- commandArgs(TRUE)
peakcaller <- args[1]
thresh <- args[2]
exp.summary <- args[3]
fig <- args[4]

cluster <- TRUE

format.annotations <- function(peakcaller,thresh,exp.summary,fig){

  peaks.overlap <- read.csv(paste0('fig1/peak_overlaps_summary_',peakcaller,'_',thresh,'.txt'),sep=' ')
  peaks.overlap$cellline1 <- gsub("_"," ",peaks.overlap$cellline1)
  peaks.overlap$cellline2 <- gsub("_"," ",peaks.overlap$cellline2)
  overlap.annotations <- read.csv(exp.summary,sep=' ')
  overlap.annotations$cell_line <- gsub("_"," ",overlap.annotations$cell_line)
  overlap.annotations$label <- gsub("_"," ",overlap.annotations$label)
  if (fig == 'xiao'){
    overlap.annotations <- overlap.annotations[which(overlap.annotations$fig == "xiao"),]
    overlap.annotations$explabel <- gsub("_"," ",paste0(overlap.annotations$label))
  }else{
    overlap.annotations <- overlap.annotations[which(overlap.annotations$fig == "1c" | overlap.annotations$fig == "1d"),]
    overlap.annotations$explabel <- gsub("_"," ",paste0(overlap.annotations$label," (",overlap.annotations$cell_line,")"))
  }
  
  rownames(overlap.annotations) <- paste(overlap.annotations$label,overlap.annotations$cell_line)
  write.table(rownames(overlap.annotations))
  peaks.overlap$explabel1 <- overlap.annotations[paste(peaks.overlap$exp1,peaks.overlap$cellline1),"explabel"]
  peaks.overlap$explabel2 <- overlap.annotations[paste(peaks.overlap$exp2,peaks.overlap$cellline2),"explabel"]

  peaks.overlap$exp1 <- as.character(peaks.overlap$exp1)
  peaks.overlap$exp2 <- as.character(peaks.overlap$exp2)
  x <- list()
  x$peaks.overlap <- peaks.overlap
  x$overlap.annotations <- overlap.annotations
  return(x)
}

#plot Fig 1c
plot.fig1c <- function(peaks.overlap,overlap.annotations,finame){
  for (cellline in unique(peaks.overlap$cellline1)){
    cell.peak.ol <- peaks.overlap[which(peaks.overlap$cellline1 == cellline & peaks.overlap$cellline2 == cellline),]
    if (dim(cell.peak.ol)[1] > 0){
      perc.overlap.mat <- matrix(data = NA,nrow=length(unique(cell.peak.ol$exp1)),ncol=length(unique(cell.peak.ol$exp2)))
      perc.overlap.mat <- as.data.frame(perc.overlap.mat)
      rownames(perc.overlap.mat) <- unique(cell.peak.ol$exp1)
      colnames(perc.overlap.mat) <- unique(cell.peak.ol$exp1)
      for (exp1 in cell.peak.ol$exp1){
        perc.overlap.mat[exp1,exp1] <- 100
        for (exp2 in cell.peak.ol[which(cell.peak.ol$exp1 == exp1),'exp2']){
          peak.overlap <- cell.peak.ol[which(cell.peak.ol$exp1 == exp1 & cell.peak.ol$exp2 == exp2),'peak_overlap'][1]
          total.peaks <- cell.peak.ol[which(cell.peak.ol$exp1 == exp1 & cell.peak.ol$exp2 == exp2),'total_peaks'][1]
          perc.overlap.mat[exp1,exp2] <- peak.overlap*100/total.peaks
        }
      }
    
      diag(perc.overlap.mat)=NA
      perc.overlap.sum <- as.vector(as.matrix(perc.overlap.mat))
      summary <- paste(summary(perc.overlap.sum))
      writeLines(paste(cellline,'1st quartile = ',summary[2],', median = ',summary[3],', 3rd quartile = ',summary[5]))
      writeLines(paste(cellline,'min = ',summary[1],', max = ',summary[6]))
    
      diag(perc.overlap.mat) = 100
      #could add the total #peaks per comparison
      pdf(paste0("fig1c_",cellline,"_peak_heatmap_",finame,".pdf"),width=3.5,height=3)
      ch = Heatmap(perc.overlap.mat,col=colorRamp2(c(0,100),c("white","#89043d")), cluster_rows = FALSE, cluster_columns = FALSE,
              column_names_side = "top",row_names_side = "left",name="% peaks\noverlapped",rect_gp = gpar(col = "white", lwd = 2),
              column_title = "Experiment 2",row_title = "Experiment 1",column_title_gp = gpar(fontsize=12),row_title_gp = gpar(fontsize=12))
      draw(ch, column_title=cellline,column_title_gp = gpar(fontsize=14))
      dev.off()
    }
  }
}

#plot Fig 1d
plot.fig1d <- function(peaks.overlap,overlap.annotations,finame,feature1="cell line",feature2="study",fig='1d'){
  for (species in c("hg38")){
    cell.lines <- overlap.annotations[which(overlap.annotations$species == species & overlap.annotations$fig == fig),"cell_line"]
    cell.peak.ol <- peaks.overlap[which(peaks.overlap$cellline1 %in% cell.lines),]
    perc.overlap.mat <- matrix(data = NA,nrow=length(unique(cell.peak.ol$explabel1)),ncol=length(unique(cell.peak.ol$explabel2)))
    perc.overlap.mat <- as.data.frame(perc.overlap.mat)
    rownames(perc.overlap.mat) <- unique(cell.peak.ol$explabel1) 
    colnames(perc.overlap.mat) <- unique(cell.peak.ol$explabel2) 
    for (explabel1 in cell.peak.ol$explabel1){
      perc.overlap.mat[explabel1,explabel1] <- 100
      for (explabel2 in cell.peak.ol[which(cell.peak.ol$explabel1 == explabel1),'explabel2']){
        peak.overlap <- cell.peak.ol[which(cell.peak.ol$explabel1 == explabel1 & cell.peak.ol$explabel2 == explabel2),'peak_overlap'][1]
        total.peaks <- cell.peak.ol[which(cell.peak.ol$explabel1 == explabel1 & cell.peak.ol$explabel2 == explabel2),'total_peaks'][1]
        perc.overlap.mat[explabel1,explabel2] <- peak.overlap*100/total.peaks
      }
    }
  }
  
  row.order <- sort(overlap.annotations[which(overlap.annotations$explabel %in% rownames(perc.overlap.mat)),]$explabel)
  perc.overlap.mat <- perc.overlap.mat[row.order,row.order]
  o.a <- overlap.annotations[which(overlap.annotations$explabel %in% row.order),]
  rownames(o.a) <- o.a$explabel
  o.a <- o.a[row.order,]
  cl <- sort(unique(gsub("_"," ",o.a$cell_line)))
  cell.colours <- list(feature1= setNames(ggpubr::get_palette(c("#3f1a1c","#d78521","#e6af2e","#ede5a6"),length(cl)),cl))
  st <- unique(o.a$label)
  if (feature2 == "gestation time"){f2.pal <- ggpubr::get_palette("Blues",length(st))
  }else{ f2.pal <- ggpubr::get_palette("lancet",length(st))}
  study.colours <- list(feature2= setNames(f2.pal,st)) 
  
  celllines = rowAnnotation(df = data.frame(feature1 = sub("_"," ",o.a$cell_line)),
                                col = cell.colours,annotation_legend_param = list(feature1 = list(title = feature1)),
                            show_annotation_name = FALSE) 
  studies = rowAnnotation(df = data.frame(feature2 = o.a$label),col = study.colours,
                          annotation_legend_param = list(feature2 = list(title = feature2)),show_annotation_name = FALSE) 
  
  pdf(paste0("fig1d_clustered_peak_heatmap_",finame,".pdf"),width=8,height=6)
  ch = Heatmap(perc.overlap.mat,col=colorRamp2(c(0,100),c("white","#89043d")), cluster_rows = TRUE, cluster_columns = TRUE, 
               column_names_side = "top",row_names_side = "left",name="% peaks\noverlapped",rect_gp = gpar(col = "white", lwd = 2),
               column_title = "Experiment 2",row_title = "Experiment 1",column_title_gp = gpar(fontsize=12),row_title_gp = gpar(fontsize=12)) 
  draw(ch + celllines + studies, heatmap_legend_side = "right")
  dev.off()
  
  perc.overlap.sum <- as.vector(as.matrix(perc.overlap.mat))
  perc.overlap.sum <- perc.overlap.sum[which(perc.overlap.sum != 100)]
  summary(perc.overlap.sum)
}

x <- format.annotations(peakcaller,thresh,exp.summary,fig)
if (fig == "xiao"){
  plot.fig1d(x$peaks.overlap,x$overlap.annotations,paste0(peakcaller,thresh),"tissue","gestation time",fig)
}else{
  #plot.fig1c(x$peaks.overlap,x$overlap.annotations,peakcaller)
  plot.fig1d(x$peaks.overlap,x$overlap.annotations,paste0(peakcaller,thresh))
}
