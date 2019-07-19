#devtools::install_github("compgenomics/MeTPeak",build_opts = c("--no-resave-data","--no-manual"))
#devtools::install_github("al-mcintyre/DEQ")
library(deq)
library(rtracklayer)
library(exomePeak)
library(MeTPeak)
library(MeTDiff)
library(UpSetR)
library(ggplot2)
library(ggpubr)
library(BSgenome)
library(doBy)

#change to wherever your mouse and human gtfs are kept
mm10.gtf <- ""
hg38.gtf <- ""

control.summary <- read.csv('summary_files/fig2_experiment_summary.txt',sep=' ')
control.summary <- control.summary[which(control.summary$control_type == "negative"),]
control.summary$studycont <- paste0(control.summary$study,"_",control.summary$control_type)

load.fi <- function(fi, anno, annotation.order,bsg){
  peaks <- deq:::import.peaks(fi,anno,annotation.order)
  peaks$peaks$seqs <- BSgenome::getSeq(bsg,peaks$peaks)
  peaks$peaks$length <- width(peaks$peaks)
  peaks$peaks$dracn <- count.matches("[AGT][AG]AC",peaks$peaks$seqs) 
  return(peaks)
}

count.matches <- function(pat, vec) sapply(regmatches(vec, gregexpr(pat, vec,perl=TRUE)), length)

for (expcont in control.summary$studycont){
  exp <- control.summary[which(control.summary$studycont == expcont),'study']
  ctype <- control.summary[which(control.summary$studycont == expcont),'control_type']
  col.data <- read.csv(paste0('fig2/samplesheet_',ctype,'_',exp,'.csv'))
  genomeTag <- as.character(control.summary[which(control.summary$study == exp & control.summary$control_type == ctype),'species'])
  if (genomeTag == "mm10"){
    gtf <- mm10.gtf
  }else {
    gtf <- hg38.gtf
  }

  input.bams <- paste0("fig2/",col.data$bam[which(col.data$IP == "input")])
  ip.bams <- paste0("fig2/",col.data$bam[which(col.data$IP == "IP")])
  input.col <- col.data[which(col.data$IP == "input"),]
  ip.col <- col.data[which(col.data$IP == "IP"),]
  metpeak.peaks <- metpeak(GENE_ANNO_GTF = gtf,IP_BAM = ip.bams,INPUT_BAM = input.bams,OUTPUT_DIR = paste0("data/",exp,"_MeTPeak_output"))
  exomepeak.peaks <- exomepeak(IP_BAM = ip.bams,INPUT_BAM = input.bams,GENE_ANNO_GTF = gtf, OUTPUT_DIR = paste0("data/",exp,"_exomePeak_output"))
  metdiff.peaks <- metdiff(GENE_ANNO_GTF = gtf, IP_BAM = ip.bams[which(ip.col$Condition == ip.col$Condition[2])],INPUT_BAM = input.bams[which(ip.col$Condition == ip.col$Condition[2])],
                           TREATED_IP_BAM = ip.bams[which(ip.col$Condition == ip.col$Condition[1])],
                           TREATED_INPUT_BAM = input.bams[which(input.col$Condition == col.data$Condition[1])],
                           OUTPUT_DIR = paste0("data/",exp,"_MeTDiff_output1")) 
}

annotation.order <- c("utr3","utr5","exon","intron","utr3*","utr5*")
peak.callers <- c("MACS2","MeTPeak","MeTDiff","exomePeak")
for (expcont in control.summary$studycont[2]){
  exp <- control.summary$study[which(control.summary$studycont == expcont)]
  label <- control.summary$control[which(control.summary$studycont == expcont)]
  col.data <- read.csv(paste0('fig2/samplesheet_negative_',exp,'.csv'))
  genomeTag <- as.character(control.summary[which(control.summary$study == exp & control.summary$control_type == "negative"),'species'])
  if (genomeTag == "mm10"){
    gtf <- mm10.gtf
    namespace <- loadNamespace("BSGenome.Mmusculus.gencode.M19")
    bsgenome <- namespace[["BSGenome.Mmusculus.gencode.M19"]]
  }else {
    gtf <- hg38.gtf
    namespace <- loadNamespace("BSgenome.Hsapiens.UCSC.hg38")
    bsgenome <- namespace[["BSgenome.Hsapiens.UCSC.hg38"]]
  }
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf,format='gtf')
  gtf.in <- rtracklayer::import(gtf)
  genenames <- unique(as.data.frame(gtf.in)[,c("gene_id","gene_name","strand")])
  anno <- list(txdb=txdb,genenames=genenames)

  macs2.file <- paste0("data/",exp,"/macs2_results/",label,"_min1.bed") 
  macs2.peaks <- load.fi(macs2.file,anno,annotation.order,bsgenome)
  metpeak.file <- paste0("data/",exp,"_MeTPeak_output/MeTPeak_output/peak.bed")
  metpeak.peaks <- load.fi(metpeak.file,anno,annotation.order,bsgenome)
  metdiff.file <- paste0("data/",exp,"_MeTDiff_output/MeTDiff_output/peak.bed")
  metdiff.peaks <- load.fi(metdiff.file,anno,annotation.order,bsgenome)
  exomepeak.file <- paste0("data/",exp,"_exomePeak_output/exomePeak_output/peak.bed")
  exomepeak.peaks <- load.fi(exomepeak.file,anno,annotation.order,bsgenome)
  
  peak.list <- c(macs2.peaks$peaks,metpeak.peaks$peaks,metdiff.peaks$peaks,exomepeak.peaks$peaks)
  all.peaks <- reduce(peak.list)
  
  #control.peaks <- GRanges(paste0(as.character(seqnames(all.peaks)),":",end(all.peaks) + 1,"-",end(all.peaks)+width(all.peaks)))
  #control.peaks <- control.peaks[which(countOverlaps(control.peaks,all.peaks) == 0),]
  #control.peaks$seqs <- BSgenome::getSeq(bsgenome,control.peaks)
  #control.peaks$length <- width(control.peaks)
  #control.peaks$dracn <- count.matches("[AGT][AG]AC",control.peaks$seqs) 
  
  expression.matrix <- as.data.frame(matrix(0,nrow=length(all.peaks),ncol=length(peak.callers)))
  colnames(expression.matrix) <- peak.callers
  expression.matrix$MACS2[which(countOverlaps(all.peaks,macs2.peaks$peaks) > 0)] <- 1 
  expression.matrix$MeTDiff[which(countOverlaps(all.peaks,metdiff.peaks$peaks) > 0)] <- 1 
  expression.matrix$MeTPeak[which(countOverlaps(all.peaks,metpeak.peaks$peaks) > 0)] <- 1 
  expression.matrix$exomePeak[which(countOverlaps(all.peaks,exomepeak.peaks$peaks) > 0)] <- 1 

  pdf(paste0(expcont,"_peakcaller_upset.pdf"),width=6,height=3)
  upset(expression.matrix,order.by="freq")
  dev.off()
  
  m2.peaks <- macs2.peaks$peaks #[which(macs2.peaks$peaks$length < 500),]
  md.peaks <- metdiff.peaks$peaks #[which(metdiff.peaks$peaks$length < 500),]
  mp.peaks <- metpeak.peaks$peaks #[which(metpeak.peaks$peaks$length < 500),]
  ep.peaks <- exomepeak.peaks$peaks #[which(exomepeak.peaks$peaks$length < 500),]
  dracn.df <- as.data.frame(rbind(cbind(m2.peaks$dracn/m2.peaks$length,"MACS2",m2.peaks$length),
                    cbind(md.peaks$dracn/md.peaks$length,"MeTDiff",md.peaks$length),
                    cbind(mp.peaks$dracn/mp.peaks$length,"MeTPeak",mp.peaks$length),
                    cbind(ep.peaks$dracn/ep.peaks$length,"exomePeak",ep.peaks$length)))
                    #(control.peaks$dracn/control.peaks$length,"control",control.peaks$length)))
  colnames(dracn.df) <- c("num.motifs","peak.caller","length")
  dracn.df$num.motifs <- as.numeric(as.character(dracn.df$num.motifs))
  summaryBy(num.motifs~peak.caller,data=dracn.df,FUN=mean)
  
  pdf(paste0("sfig1_",expcont,"_peakcaller_dracn.pdf"),width=6,height=2)
  ggplot(data=dracn.df, aes(x=peak.caller,y=num.motifs)) + geom_violin() + geom_boxplot(width=0.1) +
     xlab("Peak caller") + ylab("# DRACN motifs/\npeak length") +
     theme(panel.background = element_blank(),axis.text.x = element_text(angle=90, hjust=1)) +
     theme_classic() + theme(legend.position="none") + 
    stat_summary(fun.data = function(x) data.frame(y=0.21, label = paste("Mean=",round(mean(x),3))), geom="text")
  dev.off()
  
}
