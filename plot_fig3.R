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

source('plot_fig3.R')

args <- commandArgs(TRUE)
peakcaller <- args[1]
thresh <- args[2]
exp.summary <- args[3]

experiment.summaries <-c('kshv','hiv','hcmv','fig3')
if (peakcaller == "metdiff"){outdir="deq_output_metdiff_peaks/"}else{outdir="deq_output/"}

#set up output data frame
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
for (summary in experiment.summaries[4]){
  exp.summary <- read.csv(paste0('summary_files/',summary,'_experiment_summary.txt'),sep=' ')
  for (exp in exp.summary$study[13]){ 
    writeLines(exp)
    coldata <- read.csv(paste0('fig3/samplesheet_',exp,'.csv'))
  
    gtf <- get.gtf(exp.summary[which(exp.summary$study == exp),"species"])
    readlen <- exp.summary[which(exp.summary$study == exp),"readlen"]
    fraglen <- exp.summary[which(exp.summary$study == exp),"fraglen"] 
    pe <- exp.summary[which(exp.summary$study == exp),"pe"]
    
    
    control <- as.character(exp.summary[which(exp.summary$study == exp),"control"])
    treatment <- as.character(exp.summary[which(exp.summary$study == exp),"stimulus"])
    if (peakcaller == "metdiff"){
      peak.files <- list.files(paste0('data/',exp,'/',control,'_',treatment,'_MeTDiff_output/MeTDiff_output/'),pattern='peak.sorted.bed',full.names=TRUE)
    }else {
      peak.files <- as.vector(unique(coldata$Peaks))
    }
    outfi <- paste0(outdir,exp,"_deq_results.txt")
    
    control <- unique(coldata$Condition)[1]
    treatment <- unique(coldata$Condition)[2]
    input.bams <- paste0(as.vector(coldata$bam[which(coldata$Condition == control & coldata$IP == "input")]))
    ip.bams <- paste0(as.vector(coldata$bam[which(coldata$Condition == control & coldata$IP == "IP")]))
    treatment.input.bams <- paste0(as.vector(coldata$bam[which(coldata$Condition == treatment & coldata$IP == "input")]))
    treatment.ip.bams <- paste0(as.vector(coldata$bam[which(coldata$Condition == treatment & coldata$IP == "IP")]))
    
    results <- deq(input.bams,ip.bams,treatment.input.bams,treatment.ip.bams,peak.files,gtf,tool='deq',paired.end=pe,outfi=outfi,readlen=readlen,fraglen=fraglen)
  }
}

#evaluate tool results and plot Figure 3a
min.count <- 0
experiment.summary <- read.csv('summary_files/fig3_experiment_summary.txt',sep=' ')
for (exp in experiment.summary$study[1:13]){
  peaks <- read.csv(paste0(outdir,exp,"_deq_results.txt"),sep="\t")
  counts <- as.matrix(read.csv(paste0(outdir,exp,"_deq_results.counts.txt"),sep="\t"))

  PmG <- peaks$diff.l2fc
  s.peaks <- c()
  s.peaks.intersect <- 1:dim(peaks)[1]
    
  for (tool in tools){
    p <- peaks[,paste0(tool,'.p')]
    padj <- peaks[,paste0(tool,'.padj')]
    
    df.inds <- which(tool.results$tool == tool & tool.results$experiment == exp)
    tool.results[df.inds,"n_peaks"] <- length(p)
    tool.results[df.inds,"n_peaks_sig"] <- length(which(padj < 0.05))
    
    writeLines(paste(exp,tool,length(which(p<0.05) )*100/length(p),length(which(padj<0.05))))
    sig.peaks.this.tool <- which(padj < 0.05)
    with.min <- which(!(peaks$annot %in% c("intron","intergenic")) & abs(PmG) >= 1 & rowMin(counts) >= min.count)
    s.peaks <- union(s.peaks,intersect(sig.peaks.this.tool,with.min))
    s.peaks.intersect <- intersect(s.peaks.intersect,intersect(sig.peaks.this.tool,with.min))
    if (length(sig.peaks.this.tool) > 0){
      lowest.p <- which(padj == min(padj[sig.peaks.this.tool]))
      write.table(cbind(as.data.frame(peaks[lowest.p,])[,c("seqnames","start","end","main.gene","width","strand","annot","diff.l2fc")],padj[lowest.p]),paste0(exp,"_",tool,"_sig_peaks.bed"),quote=FALSE,sep="\t",row.names = FALSE, col.names = FALSE)
      write.table(unique(as.data.frame(peaks[sig.peaks.this.tool,])$maingene),paste0(exp,"_",tool,"_sig_genes.bed"),quote=FALSE,sep="\t",row.names = FALSE, col.names = FALSE)
    }
  }
  tool <- "intersect"
  df.inds <- which(tool.results$tool == tool & tool.results$experiment == exp)
  tool.results[df.inds,"n_peaks_sig"] <- length(s.peaks.intersect)
  tool <- "filtered"
  df.inds <- which(tool.results$tool == tool & tool.results$experiment == exp)
  tool.results[df.inds,"n_peaks_sig"] <- length(s.peaks)
  tool <- "original"
  df.inds <- which(tool.results$tool == tool & tool.results$experiment == exp)
  tool.results[df.inds,"n_peaks_sig"] <- experiment.summary$reported_peaks[which(experiment.summary$study == exp)]
  if (length(s.peaks) > 0){
    write.table(as.data.frame(peaks[s.peaks,])[,c("seqnames","start","end","main.gene","width","strand","annot","diff.l2fc")],paste0(exp,"_",peakcaller,"_sig_peaks.bed"),quote=FALSE,sep="\t",row.names = FALSE, col.names = FALSE)
    write.table(unique(as.data.frame(peaks[s.peaks,])$main.gene),paste0(exp,"_",peakcaller,"_sig_genes.bed"),quote=FALSE,sep="\t",row.names = FALSE, col.names = FALSE)
    
  }
}

summary <- tool.results[complete.cases(tool.results[,c('tool','experiment','n_peaks_sig','n_reps','Condition')]),]
filtered.summary <- summary[which(summary$tool == "filtered"),]
condition.order <- filtered.summary[order(filtered.summary$n_peaks_sig, -filtered.summary$n_reps),'Condition']
summary$Condition <- factor(summary$Condition,levels = condition.order)
summary$tool <- factor(summary$tool,levels = more.tools)
summary$n_peaks_sig <- summary$n_peaks_sig + 1

pdf(paste0("fig3a_number_sig_peaks_min",min.count,"_",peakcaller,".pdf"),width = 6, height=4)
ggplot(data=summary, aes(x=Condition, y=n_peaks_sig, fill=tool)) +
  geom_bar(stat="identity", position=position_dodge()) + ylab("# peaks changed + 1 (log10 scale)") + 
  facet_grid(~ Condition, space="free_x", scales="free_x", switch="x") + 
  xlab("") + scale_y_continuous(trans='log10') +
  scale_fill_manual(values=tool.cols) + theme(panel.grid.major.y = element_line(colour = "lightgrey"), panel.grid.major.x = element_blank(),
                                              panel.background = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.background = element_blank(),strip.text.x = element_blank(),panel.spacing = unit(0.2, "lines"))
dev.off()


#reanalyze experiments under similar conditions (HIV, KSHV, dsDNA)

reanalyses <- c("hiv","kshv","hcmv")
for (condition in reanalyses[1]){
  experiment.summary <- read.csv(paste0('summary_files/',condition,'_experiment_summary.txt'),sep=' ')
  sig.genes <- list()
  sig.peaks <- list()
  for (exp in experiment.summary$study){
    peaks <- read.csv(paste0(outdir,exp,"_deq_results.txt"),sep="\t")
    PmG <- peaks$diff.l2fc
    genes <- c()
    s.peaks <- c()
    for (tool in tools){
      p <- peaks[,paste0(tool,'.p')]
      padj <- peaks[,paste0(tool,'.padj')]
      writeLines(paste(exp,tool,length(which(p<0.05) )*100/length(p)))
      genes <- union(genes,peaks$main.gene[which(padj < 0.05 & abs(PmG) >= 1)])
      s.peaks <- union(s.peaks,which(padj < 0.05 & abs(PmG) >= 1))
    }
    sig.genes[[exp]] <- genes
    sig.peaks[[exp]] <- GenomicRanges::makeGRangesFromDataFrame(peaks[s.peaks,],keep.extra.columns = TRUE)
    write.table(as.data.frame(sig.peaks[[exp]])[,c("seqnames","start","end","main.gene","width","strand","annot","diff.l2fc")],paste0(exp,"_sig_peaks.bed"),quote=FALSE,sep="\t",row.names = FALSE, col.names = FALSE)
  }
  if (condition == "hcmv"){
    hcmv.intersect <- GenomicRanges::intersect(sig.peaks[["winkler_hcmv_2018"]],
                                sig.peaks[["rubio_hcmv_2018"]])
    write.table(as.data.frame(hcmv.intersect)[,c("seqnames","start","end")],paste0("hcmv_sig_peaks.bed"),quote=FALSE,sep="\t",row.names = FALSE, col.names = FALSE)
  }
  if (condition == "hiv"){
    hiv.intersect <- GenomicRanges::intersect(sig.peaks[["tirumuru_hiv_2016"]],
                                               sig.peaks[["lichinchi_hiv_2016"]])
    write.table(as.data.frame(hiv.intersect)[,c("seqnames","start","end")],paste0("hcmv_sig_peaks.bed"),quote=FALSE,sep="\t",row.names = FALSE, col.names = FALSE)
  }
  if (condition == "kshv"){
    kshv.intersect.1 <- GenomicRanges::intersect(sig.peaks[["hesser_kshv_2018"]],sig.peaks[["tan_kshv_2018"]])
    if (length(kshv.intersect.1) > 0){
      write.table(as.data.frame(kshv.intersect.1)[,c("seqnames","start","end")],paste0("kshv_sig_peaks.1.bed"),quote=FALSE,sep="\t",row.names = FALSE, col.names = FALSE)
    }
    kshv.intersect.2 <- GenomicRanges::intersect(sig.peaks[["hesser_kshv_2018"]],sig.peaks[["tan_kshv_2018a"]])
    if (length(kshv.intersect.2) > 0){
      write.table(as.data.frame(kshv.intersect.2)[,c("seqnames","start","end")],paste0("kshv_sig_peaks.2.bed"),quote=FALSE,sep="\t",row.names = FALSE, col.names = FALSE)
    }
    kshv.intersect.3 <- GenomicRanges::intersect(sig.peaks[["hesser_kshv_2018"]],sig.peaks[["tan_kshv_2018b"]])
    if (length(kshv.intersect.3) > 0){
      write.table(as.data.frame(kshv.intersect.3)[,c("seqnames","start","end")],paste0("kshv_sig_peaks.3.bed"),quote=FALSE,sep="\t",row.names = FALSE, col.names = FALSE)
    }
  }
}

#for testing
# txdb <- GenomicFeatures::makeTxDbFromGFF(gtf,format='gtf')
# gtf.in <- rtracklayer::import(gtf)
# genenames <- unique(as.data.frame(gtf.in)[,c("gene_id","gene_name","strand")])
# anno <- list(txdb=txdb,genenames=genenames)
# annotation.order <- c("utr3","utr5","exon","intron","utr3*","utr5*")
# 
# peaks <- deq:::import.peaks(peak.files,anno,annotation.order)
# 
# all.bams <- c(input.bams,ip.bams,treatment.input.bams,treatment.ip.bams)
# peaks <- deq:::count.reads(peaks,all.bams,pe,fraglen-readlen)
