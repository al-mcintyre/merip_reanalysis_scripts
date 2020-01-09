library(S4Vectors)
library(dplyr)
library(deq)
library(doBy)
library(moments)
library(scatterplot3d)
library(ggplot2)
library(ggpubr) #for ggscatter
library(ggfortify) #needed for autoplot
library(Biobase) #for rowMin

source('plot_fig3_helper.R')

args <- commandArgs(TRUE)
peakcaller <- args[1] #macs
thresh <- as.numeric(args[2])
exp.summary <- args[3] #exp_summary_mergedpeaks.txt
mm10.gtf <- args[4]
hg38.gtf <- args[5]
run.deq <- args[6]

experiment.summaries <-c('kshv','hiv','dsDNA')
outdir <- paste0('fig4/deq_output_mergedpeaks/')

#set up output data frame
experiment.summary <- read.csv(exp.summary,sep=' ')
nexps <- dim(experiment.summary)[1]
tools <- c("deseq2","edger","qnb")
more.tools <- c("original",tools,"intersect","filtered")
output.columns <- c("tool","experiment","n_peaks","n_peaks_sig","percent_under0.05","n_reps")
tool.results <- as.data.frame(matrix(ncol=length(output.columns), nrow=nexps*length(more.tools)))
colnames(tool.results) <- output.columns
tool.results$tool <- sort(rep(more.tools,nexps))
tool.results$experiment <- rep(experiment.summary$study,length(more.tools))
tool.results$n_reps <- rep(experiment.summary$reps,length(more.tools))
tool.results$Condition <- rep(experiment.summary$condition,length(more.tools))

#plot colours
tool.cols <- c("#707070","#66D7D1","#1b9aaa","#db7b78","#a9a587","#721652")
tool.colours <- data.frame(colour=tool.cols)
row.names(tool.colours) <- more.tools

#run DESeq2, edgeR, and QNB 
if (run.deq == "TRUE"){
for (summary in experiment.summaries){
  exp.summary <- experiment.summary[which(experiment.summary$fig == summary),] #%in% experiment.summaries),]
  for (exp in exp.summary$study){ 
    writeLines(exp)
  
    subexpsum <- experiment.summary[which(experiment.summary$study == exp),]
    exp <- subexpsum$study

    gtf <- get.gtf(subexpsum$species,mm10.gtf,hg38.gtf)
    readlen <- subexpsum$readlen
    fraglen <- subexpsum$fraglen
    pe <- subexpsum$pe
    outfi <- paste0(outdir,exp,"_deq_results.txt")

    control <- as.character(subexpsum$control)
    stimulus <- as.character(subexpsum$stimulus)

    input.bams <- list.files(paste0(exp,'/alignments'),pattern=paste0(control,'_input_[1-9].star.sorted.bam$'),full.names=TRUE)
    ip.bams <- list.files(paste0(exp,'/alignments'),pattern=paste0(control,'_IP_[1-9].star.sorted.bam$'),full.names=TRUE)
    stimulus.input.bams <- list.files(paste0(exp,'/alignments'),pattern=paste0(stimulus,'_input_[1-9].star.sorted.bam$'),full.names=TRUE)
    stimulus.ip.bams <- list.files(paste0(exp,'/alignments'),pattern=paste0(stimulus,'_IP_[1-9].star.sorted.bam$'),full.names=TRUE)

    if (peakcaller == "merged"){
        peak.files <- c(paste0('merged_',summary,'_peaks.bed'))
    }

    results <- deq(input.bams,ip.bams,stimulus.input.bams,stimulus.ip.bams,peak.files,gtf,tool='deq',paired.end=pe,outfi=outfi,readlen=readlen,fraglen=fraglen)
  }
}
}

for (summary in experiment.summaries){
  exp.summary <- experiment.summary[which(experiment.summary$fig == summary),]
  diffs <- list()
  counts.over.min <- list()
  exp1 <- paste(exp.summary$study[1])
  label1 <- paste(exp.summary$label[1])
  exp2 <- paste(exp.summary$study[2])
  label2 <- paste(exp.summary$label[2])
  writeLines(paste(exp1,exp2,label1,label2))
  for (exp in exp.summary$study){
    peaks <- read.csv(paste0(outdir,exp,"_deq_results.txt"),sep="\t")
    counts <- as.matrix(read.csv(paste0(outdir,exp,"_deq_results.counts.txt"),sep="\t"))
    #write.table(head(rowMin(counts) >= thresh))
    #write.table( grepl("input", colnames(counts)))
    #write.table(head(counts))
    counts.over.min[[exp]] <- rowMin(counts) >= thresh # rowMedians(counts[, grepl("input", colnames(counts))]) >= thresh
    diffs[[exp]] <- peaks$diff.l2fc
  }
  
  diffs.df <- data.frame(diffs)
  diffs.df2 <- data.frame(apply(diffs.df, 2, function(x) as.numeric(as.character(x))))
  counts.min.df <- data.frame(counts.over.min)
  #write.table(head(counts.min.df))
  diffs.df2 <- diffs.df2[which(counts.min.df[,exp1] & counts.min.df[,exp2]),]
  outfi.df <- cbind(peaks[which(counts.min.df[,exp1] & counts.min.df[,exp2]),c("main.gene","annot")], diffs.df2)
  write.table(outfi.df,paste0('fig4/',summary,'_min',thresh,'.txt'),quote=FALSE,row.names=FALSE,col.names=TRUE)
  #write.table(head(diffs.df2,10))
  #writeLines(paste(exp2))
  #writeLines(class(diffs.df2[,paste(exp2)]))

  pdf(paste0('fig4/sfig4_',summary,'_scatter_min',thresh,'.pdf'),height=4,width=4)
  g <- ggscatter(diffs.df2,x=exp1,y=exp2, add="reg.line", add.params = list(color = "black"),
                      conf.int = FALSE,cor.coef = TRUE, cor.method = "pearson",
                      xlab=paste0("Peak - gene L2FC, ",label1), ylab=paste0("Peak - gene L2FC, ",label2), color="#89043d",
                      alpha=0.4,legend.title = '',legend='top',font.legend = c(8, "plain", "black"))
  print(g)
  dev.off()

}
