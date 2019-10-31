library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(reshape2)
library(dplyr)

#plot SFig 1e number of reads vs number of peaks

plot.sfig1e <- function(peakcaller,read.label){
    finame <- paste0('fig1/reads_vs_',peakcaller,'_peaks_',read.label,'.txt')
    reads.v.peaks <- read.csv(finame,sep=" ")
    reads.v.peaks$label <- gsub('_',' ',reads.v.peaks$label)

    pdf(paste0('fig1/sfig1d_',read.label,'_',peakcaller,'_reads_vs_peaks.pdf'),height=6,width=7)
    g <- ggscatter(reads.v.peaks,x="reads",y="peaks", add="reg.line",add.params = list(color = "black"),
          conf.int = FALSE, cor.coef = TRUE, cor.method = "pearson",
          xlab=paste("#",read.label,"reads"), ylab="# peaks",color="label",shape="label",alpha=1,
          legend.title = '',legend='top',font.legend = c(7, "plain", "black")) + scale_shape_manual(values=rep(c(8,seq(15,18)),8))
    grid.draw(g)
    dev.off()  
}

args <- commandArgs(TRUE)
peakcaller <- args[1]
read.label <- args[2]
plot.sfig1e(peakcaller,read.label)
