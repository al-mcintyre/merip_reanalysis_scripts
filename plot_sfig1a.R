#requires MACS2 and MeTDiff results for negative control data sets -- runs MeTPeak and exomePeak
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

args <- commandArgs(TRUE)
mm10.gtf <- args[1]
hg38.gtf <- args[2]
exp.summary <- args[3]
run.peakcallers <- TRUE

source('rhelper.R')

control.summary <- read.csv(exp.summary,sep=' ')
control.summary <- control.summary[which(control.summary$control_type == "negative"),]
control.summary$studycont <- paste0(control.summary$study,"_",control.summary$control_type)

count.matches <- function(pat, vec) sapply(regmatches(vec, gregexpr(pat, vec,perl=TRUE)), length)

if (run.peakcallers == "TRUE"){
for (expcont in control.summary$studycont){
  exp <- control.summary[which(control.summary$studycont == expcont),'study']
  ctype <- control.summary[which(control.summary$studycont == expcont),'control_type']
  control <- control.summary[which(control.summary$studycont == expcont),'control']
  genomeTag <- as.character(control.summary[which(control.summary$study == exp & control.summary$control_type == ctype),'species'])
  if (genomeTag == "mm10"){
    gtf <- mm10.gtf
  }else {
    gtf <- hg38.gtf
  }

  input.bams <- list.files(paste0(exp,'/alignments'),pattern=paste0("^",control,"_input_[1-9].star.sorted.bam$"),full.names=TRUE)
  ip.bams <- list.files(paste0(exp,'/alignments'),pattern=paste0("^",control,"_IP_[1-9].star.sorted.bam$"),full.names=TRUE)
  metpeak.peaks <- metpeak(GENE_ANNO_GTF = gtf,IP_BAM = ip.bams,INPUT_BAM = input.bams,OUTPUT_DIR = paste0(exp,"/MeTPeak_output"))
  exomepeak.peaks <- exomepeak(IP_BAM = ip.bams,INPUT_BAM = input.bams,GENE_ANNO_GTF = gtf, OUTPUT_DIR = paste0(exp,"/exomePeak_output"))
}
}

annotation.order <- c("utr3","utr5","exon","intron","utr3*","utr5*")
peak.callers <- c("MACS2","MeTPeak","MeTDiff","exomePeak")
for (expcont in control.summary$studycont[2]){
  exp <- control.summary$study[which(control.summary$studycont == expcont)]
  control <- control.summary$control[which(control.summary$studycont == expcont)]
  genomeTag <- as.character(control.summary[which(control.summary$study == exp & control.summary$control_type == "negative"),'species'])
  if (genomeTag == "mm10"){
    gtf <- mm10.gtf
    namespace <- loadNamespace("BSgenome.Mmusculus.UCSC.mm10")
    bsgenome <- namespace[["BSgenome.Mmusculus.UCSC.mm10"]]
    #namespace <- loadNamespace("BSGenome.Mmusculus.gencode.M19")
    #bsgenome <- namespace[["BSGenome.Mmusculus.gencode.M19"]]
  }else {
    gtf <- hg38.gtf
    namespace <- loadNamespace("BSgenome.Hsapiens.UCSC.hg38")
    bsgenome <- namespace[["BSgenome.Hsapiens.UCSC.hg38"]]
  }
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf,format='gtf')
  gtf.in <- rtracklayer::import(gtf)
  genenames <- unique(as.data.frame(gtf.in)[,c("gene_id","gene_name","strand")])
  anno <- list(txdb=txdb,genenames=genenames)

  macs2.file <- paste0(exp,"/macs2_results/",control,"_min1.bed") 
  macs2.peaks <- load.fi(macs2.file,anno,annotation.order,bsgenome)
  metpeak.file <- paste0(exp,"/MeTPeak_output/peak.bed")
  metpeak.peaks <- load.fi(metpeak.file,anno,annotation.order,bsgenome)
  metdiff.file <- paste0(exp,"/",control,"_MeTDiff_output/peak.bed")
  metdiff.peaks <- load.fi(metdiff.file,anno,annotation.order,bsgenome)
  exomepeak.file <- paste0(exp,"/exomePeak_output/peak.bed")
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

  pdf(paste0('fig1/sfig1a_',expcont,"_peakcaller_upset.pdf"),width=6,height=3)
  upset(expression.matrix,order.by="freq")
  dev.off()
  
  m2.peaks <- macs2.peaks$peaks 
  md.peaks <- metdiff.peaks$peaks
  mp.peaks <- metpeak.peaks$peaks 
  ep.peaks <- exomepeak.peaks$peaks
  dracn.df <- as.data.frame(rbind(cbind(m2.peaks$dracn/m2.peaks$length,"MACS2",m2.peaks$length),
                    cbind(md.peaks$dracn/md.peaks$length,"MeTDiff",md.peaks$length),
                    cbind(mp.peaks$dracn/mp.peaks$length,"MeTPeak",mp.peaks$length),
                    cbind(ep.peaks$dracn/ep.peaks$length,"exomePeak",ep.peaks$length)))
                    #(control.peaks$dracn/control.peaks$length,"control",control.peaks$length)))
  colnames(dracn.df) <- c("num.motifs","peak.caller","length")
  dracn.df$num.motifs <- as.numeric(as.character(dracn.df$num.motifs))
  summaryBy(num.motifs~peak.caller,data=dracn.df,FUN=mean)
  
  pdf(paste0("fig1/sfig1a_",expcont,"_peakcaller_dracn.pdf"),width=6,height=2)
  ggplot(data=dracn.df, aes(x=peak.caller,y=num.motifs)) + geom_violin() + geom_boxplot(width=0.1) +
     xlab("Peak caller") + ylab("# DRACN motifs/\npeak length") +
     theme(panel.background = element_blank(),axis.text.x = element_text(angle=90, hjust=1)) +
     theme_classic() + theme(legend.position="none") + 
    stat_summary(fun.data = function(x) data.frame(y=0.21, control = paste("Mean=",round(mean(x),3))), geom="text")
  dev.off()
  
  set.seed(741)
  unique.peaks <- list()
  if (expcont == "neg_control_mocks_negative"){
    for (col in c('MACS2','MeTPeak','MeTDiff','exomePeak')){
      u.p <- all.peaks[which(expression.matrix[col] == 1 & rowSums(expression.matrix) == 1),]
      #unique.peaks[[col]] <- as.data.frame(deq:::annotate.peaks(u.p,anno,annotation.order)$peaks)
      write.table(dplyr::sample_n(unique.peaks[[col]][which(unique.peaks[[col]]$annot %in% c('utr3','utr5','exon')),],5),
                         file=paste0('Huh7_',col,'_unique_peaks.txt'),quote = FALSE)
    }
  }
  
}
