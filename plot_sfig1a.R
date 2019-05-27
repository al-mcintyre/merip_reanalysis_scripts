#library(devtools)
#install_github("compgenomics/MeTPeak",build_opts = c("--no-resave-data","--no-manual"))
library(CLIPanalyze)
library(rtracklayer)
library(exomePeak)
library(MeTPeak)
library(MeTDiff)
library(UpSetR)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSGenome.Mmusculus.gencode.M19)
library(ggplot2)
library(ggpubr)
library(doBy)

setwd('/Volumes/Markarian_501/merip_reanalysis/fig1b/')
control.summary <- read.csv('../fig2/control_summary.txt',sep=' ')
control.summary <- control.summary[which(control.summary$control_type == "negative"),]
control.summary$studycont <- paste0(control.summary$study,"_",control.summary$control_type)
source('../fig2/plot_fig2_helper.R')

load.fi <- function(fi,cdata,bsgenome){
  if(file.size(fi) > 0){
    peaks <- import.peaks(fi)
    peaks <- finish.annotation(cdata,peaks,bsg=bsgenome)
  } else{
    peaks <- GRanges()
  }
  return(peaks)
}

for (expcont in control.summary$studycont[1]){
  exp <- control.summary[which(control.summary$studycont == expcont),'study']
  ctype <- control.summary[which(control.summary$studycont == expcont),'control_type']
  col.data <- read.csv(paste0('../fig2/samplesheet_',ctype,'_',exp,'.csv'))
  genomeTag <- as.character(control.summary[which(control.summary$study == exp & control.summary$control_type == ctype),'species'])
  if (genomeTag == "mm10"){
    gtf <- "/Volumes/Markarian_501/merip_reanalysis/fig2/gencode.vM19.annotation.gtf"
  }else {
    gtf <- "/Volumes/Markarian_501/20180228_flaviviridae/gencode.v25.annotation.gtf"
  }

  input.bams <- paste0("../fig2/",col.data$bam[which(col.data$IP == "input")])
  ip.bams <- paste0("../fig2/",col.data$bam[which(col.data$IP == "IP")])
  input.col <- col.data[which(col.data$IP == "input"),]
  ip.col <- col.data[which(col.data$IP == "IP"),]
  metpeak.peaks <- metpeak(GENE_ANNO_GTF = gtf,IP_BAM = ip.bams,INPUT_BAM = input.bams,OUTPUT_DIR = paste0(exp,"_MeTPeak_output"))
  exomepeak.peaks <- exomepeak(IP_BAM = ip.bams,INPUT_BAM = input.bams,GENE_ANNO_GTF = gtf, OUTPUT_DIR = paste0(exp,"_exomePeak_output"))
  metdiff.peaks <- metdiff(GENE_ANNO_GTF = gtf, IP_BAM = ip.bams,INPUT_BAM = input.bams,
                           TREATED_IP_BAM = ip.bams[which(ip.col$Condition == ip.col$Condition[1])],
                           TREATED_INPUT_BAM = input.bams[which(input.col$Condition == col.data$Condition[1])],
                           OUTPUT_DIR = paste0(exp,"_MeTDiff_output")) 
  
  #then take exonic peaks only for macs2 for comparison
}

annotation.order <- c("utr3","utr5","exon","intron","utr3*","utr5*")
peak.callers <- c("MACS2","MeTPeak","MeTDiff","exomePeak")
for (expcont in control.summary$studycont[1]){
  exp <- control.summary$study[which(control.summary$studycont == expcont)]
  label <- control.summary$control[which(control.summary$studycont == expcont)]
  col.data <- read.csv(paste0('../fig2/samplesheet_negative_',exp,'.csv'))
  genomeTag <- as.character(control.summary[which(control.summary$study == exp & control.summary$control_type == "negative"),'species'])
  if (genomeTag == "mm10"){
    gtf <- "/Volumes/Markarian_501/merip_reanalysis/fig2/gencode.vM19.annotation.gtf"
    namespace <- loadNamespace("BSGenome.Mmusculus.gencode.M19")
    bsgenome <- namespace[["BSGenome.Mmusculus.gencode.M19"]]
  }else {
    gtf <- "/Volumes/Markarian_501/20180228_flaviviridae/gencode.v25.annotation.gtf"
    bsgenome <- load.bsgenome(genomeTag)
  }
  annot <- loadAnnot(genomeTag)
  
  macs2.file <- paste0("macs2_peaks/",label,"_min1.bed") 
  macs2.peaks <- load.fi(macs2.file,col.data,bsgenome)
  metpeak.file <- paste0(exp,"_MeTPeak_output/MeTPeak_output/peak.bed")
  metpeak.peaks <- load.fi(metpeak.file,col.data,bsgenome)
  metdiff.file <- paste0(exp,"_MeTDiff_output/MeTDiff_output/peak.bed")
  metdiff.peaks <- load.fi(metdiff.file,col.data,bsgenome)
  exomepeak.file <- paste0(exp,"_exomePeak_output/exomePeak_output/peak.bed")
  exomepeak.peaks <- load.fi(exomepeak.file,col.data,bsgenome)
  
  peak.list <- c(macs2.peaks$peaks,metpeak.peaks$peaks,metdiff.peaks$peaks,exomepeak.peaks$peaks)
  all.peaks <- reduce(peak.list)
  
  control.peaks <- GRanges(paste0(as.character(seqnames(all.peaks)),":",end(all.peaks) + 1,"-",end(all.peaks)+width(all.peaks)))
  control.peaks <- control.peaks[which(countOverlaps(control.peaks,all.peaks) == 0),]
  control.peaks <- clip.annotate(control.peaks,annot,annotation.order)
  control.peaks <- finish.annotation(col.data,control.peaks,bsgenome)
  
  expression.matrix <- as.data.frame(matrix(0,nrow=length(all.peaks),ncol=length(peak.callers)))
  colnames(expression.matrix) <- peak.callers
  expression.matrix$MACS2[which(countOverlaps(all.peaks,macs2.peaks$peaks) > 0)] <- 1 
  expression.matrix$MeTDiff[which(countOverlaps(all.peaks,metdiff.peaks$peaks) > 0)] <- 1 
  expression.matrix$MeTPeak[which(countOverlaps(all.peaks,metpeak.peaks$peaks) > 0)] <- 1 
  expression.matrix$exomePeak[which(countOverlaps(all.peaks,exomepeak.peaks$peaks) > 0)] <- 1 

  pdf(paste0(expcont,"_peakcaller_upset.pdf"),width=6,height=3)
  upset(expression.matrix,order.by="freq")
  dev.off()
  
  dracn.df <- as.data.frame(rbind(cbind(macs2.peaks$peaks$dracn/macs2.peaks$peaks$length,"MACS2",macs2.peaks$peaks$length),
                    cbind(metdiff.peaks$peaks$dracn/metdiff.peaks$peaks$length,"MeTDiff",metdiff.peaks$peaks$length),
                    cbind(metpeak.peaks$peaks$dracn/metpeak.peaks$peaks$length,"MeTPeak",metpeak.peaks$peaks$length),
                    cbind(exomepeak.peaks$peaks$dracn/exomepeak.peaks$peaks$length,"exomePeak",exomepeak.peaks$peaks$length)))#,
                    #cbind(control.peaks$peaks$dracn/control.peaks$peaks$length,"control")))
  colnames(dracn.df) <- c("num.motifs","peak.caller","length")
  dracn.df$num.motifs <- as.numeric(as.character(dracn.df$num.motifs))
  #summaryBy(num.motifs~peak.caller,data=dracn.df,FUN=mean)
  
  pdf(paste0("sfig1_",expcont,"_peakcaller_dracn.pdf"),width=6,height=2)
  ggplot(data=dracn.df, aes(x=peak.caller,y=num.motifs)) + geom_violin() + geom_boxplot(width=0.1) +
     xlab("Peak caller") + ylab("# DRACN motifs/\npeak length") +
     theme(panel.background = element_blank(),axis.text.x = element_text(angle=90, hjust=1)) +
     theme_classic() + theme(legend.position="none") + 
    stat_summary(fun.data = function(x) data.frame(y=0.21, label = paste("Mean=",round(mean(x),3))), geom="text")
  dev.off()
  
}
