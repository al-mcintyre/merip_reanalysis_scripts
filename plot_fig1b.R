library(GenomicRanges)
library(ggplot2)
library(rtracklayer)
library(reshape2)
library(grid)

plot.fig1b <- function(dir,control){
  dataset <- strsplit(dir,"_")[[1]][1]
  files <- intersect(list.files(path=paste0(dir,"/macs2_results"),pattern=control,full.names=TRUE),
                     list.files(path=paste0(dir,"/macs2_results"),pattern = "_[1-9].macs2_peaks.narrowPeak.peaks.bed$",full.names = TRUE))
  write.table(files)
  max.reps <- length(files)
  peak.totals <- matrix(data = NA,nrow=max.reps,ncol=length(files))
  
  set.seed(174)
  
  for (n in c(1:length(files))){
    combns <- combn(files,n)
    if (n >1&n<length(files)){ 
      ncombns <- min(max.reps,dim(combns)[2])  
      len <- dim(combns)[2]
      file.combns <- combns[,sample(c(1:len),ncombns)]
    } else{
      if (n == 1){
        file.combns <- files
        ncombns <- min(max.reps,length(files))
      } else { 
        file.combns <- as.data.frame(files)
        ncombns <- 1
      }
    }
  
    for (m in c(1:ncombns)){
      if (is.null(dim(file.combns)[1])){
        select.samples <- file.combns[m]
      } else { select.samples <- file.combns[,m] }
      for (file in select.samples){
        if (exists("m2peaks")){
          m2peaks <- union(m2peaks,import(file,format="BED"))
        }else{
          m2peaks <- import(file,format="BED")
        }
      }
      m2peaks <- reduce(m2peaks)
      m2peaks <- keepStandardChromosomes(m2peaks,pruning.mode = "coarse")
      peak.totals[m,n] <- length(m2peaks)
      m2peaks <- GRanges()
    }
  }
  
  peak.tots <- as.data.frame(t(peak.totals))
  peak.tots$Replicates <- c(1:length(files))
  peak.tots <- melt(peak.tots,id="Replicates")
  peak.tots$Replicates <- factor(peak.tots$Replicates)
  max.peaks <- max(peak.tots$value[!is.na(peak.tots$value)])
  peak.tots$percent.value <- 100*peak.tots$value/max.peaks
  
  for (reps in c(1,2,3,6,7)){
    median.perc.found <- median(peak.tots[which(peak.tots$Replicates == reps & !is.na(peak.tots$value)),"value"])*100/max.peaks
    writeLines(paste('median of',median.perc.found,'% total peaks found with',reps,'replicates'))
  }
  
  g.count <- ggplot(peak.tots,aes(x=Replicates,y=percent.value)) + geom_boxplot(color="#89043d", fill="#be7695", outlier.shape=NA) + geom_jitter(position=position_jitter(width=.2, height=0),col="black") + 
    xlab("# Replicates") + ylab("% Peaks") + theme(panel.background = element_blank()) + 
    scale_y_continuous(limits = c(0,100),breaks = seq(0,100,10),sec.axis = sec_axis(~ . /100*max.peaks,name="# Peaks"))
  pdf(paste0('fig1/fig1b_npeaks_by_rep_',dataset,'.pdf'),height=3.5,width=4)
  grid.draw(g.count)
  dev.off()
}

args <- commandArgs(TRUE)
dir <- args[1]
control <- args[2]
writeLines(paste('Plotting',dir,'Figure 1b'))
plot.fig1b(dir,control)
