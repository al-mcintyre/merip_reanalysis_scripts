library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(reshape2)
library(dplyr)

args <- commandArgs(TRUE)
mm10.gtf <- args[1]
hg38.gtf <- args[2]
peakcaller <- args[3]    
exp.summary <- args[4]

source('rhelper.R')
source('plot_fig3_helper.R')

#plot SFig 1g DRAC enrichment
annotation.order <- c("utr3","utr5","exon","intron","utr3*","utr5*")
control.summary <- read.csv(exp.summary,sep=' ')
control.summary <- control.summary[which(control.summary$fig == 1 & control.summary$control_type == "negative"),]
control.summary$studycont <- paste0(control.summary$study,"_",control.summary$control_type)
output.columns <- c("experiment","replicates","percent.dracn","num.reads")

for (expcont in control.summary$studycont){
  exp <- control.summary[which(control.summary$studycont == expcont),]$study
  ctype <- control.summary[which(control.summary$studycont == expcont),]$control_type
  genomeTag <- as.character(control.summary[which(control.summary$study == exp & control.summary$control_type == ctype),'species'])
  if (genomeTag == "mm10"){
    gtf <- mm10.gtf
    #namespace <- loadNamespace("BSGenome.Mmusculus.gencode.M19")
    namespace <- loadNamespace("BSgenome.Mmusculus.UCSC.mm10")
    bsgenome <- namespace[["BSgenome.Mmusculus.UCSC.mm10"]]
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
  
  label <- control.summary[which(control.summary$studycont == expcont),'control']
  dracn.summary <- as.data.frame(matrix(ncol=length(output.columns), nrow=0))
  colnames(dracn.summary) <- output.columns
  drac.per.len <- as.data.frame(matrix(ncol=3,nrow=0))
  nreps <- control.summary[which(control.summary$studycont == expcont),'reps']
  
  for (replicate in seq(nreps)){
    macs2.file <- paste0(exp,"/macs2_results/",label,"_replicate_intersect_",replicate,".bed")
    if(file.size(macs2.file) > 0){
      peaks <- load.fi(macs2.file,anno,annotation.order,bsgenome)
      num.peaks <- length(peaks$peaks)
      num.dracn <- length(which(peaks$peaks$dracn > 0)) # == "DRACN"))
      perc <- num.dracn*100/num.peaks
      drac.per.len <- rbind(drac.per.len,cbind(peaks$peaks$dracn/peaks$peaks$length,replicate,peaks$peaks$length))
    } else{
      num.peaks <- 0
      num.dracn <- 0
      perc <- 0
    }
    if (num.peaks > 0){
      sub.results <- data.frame("experiment"=c(as.character(exp)),"replicates"=c(replicate),
                                 "percent dracn"=c(perc),"peaks"=c(num.peaks))
      dracn.summary <- rbind(dracn.summary,sub.results)
    }
  }
  nreps <- dim(dracn.summary)[1]
  
  writeLines(paste('% DRAC with 1 replicate =',dracn.summary[1,3]))
  writeLines(paste('% DRAC with 2 replicates =',dracn.summary[2,3]))
  writeLines(paste('% DRAC with',nreps,'replicates =',dracn.summary[nreps,3]))          
  dracn.summary2 <- melt(dracn.summary, id=c("experiment","replicates")) 

  colnames(drac.per.len) <- c("value","replicates","length")  
  drac.per.len$value <- as.numeric(as.character(drac.per.len$value))
  drac.per.len$replicates <- as.factor(drac.per.len$replicates)
  drac.per.len$plotType <- 1
  drac.per.len <- drac.per.len[,!is.na(colnames(drac.per.len))]
  
  dracn.summary3 <- dracn.summary[,c("peaks","replicates","percent.dracn")]
  colnames(dracn.summary3) <- c("value","replicates","length") 
  dracn.summary3$plotType <- 2
  drac.combined.summary <- rbind(drac.per.len,dracn.summary3)
  
  pdf(paste0('fig1/sfig1g_',exp,'_macs2_motif_enrichment.pdf'),height=3.2,width=4.3)
  g <- ggplot(data=drac.combined.summary, aes(x=replicates,y=value)) +
    xlab("Replicates") + ylab(NULL) + 
    facet_wrap(plotType~., ncol=1,scales = "free_y",strip.position = "left", shrink = FALSE,
               labeller = as_labeller(c("2" = "Replicate overlap", "1"= "# DRAC/peak length") )) +
    geom_violin(data = filter(drac.combined.summary,plotType==1),draw_quantiles=c(0.5),fill="#b2527b") + 
    geom_bar(data = filter(drac.combined.summary,plotType==2), stat = "identity",fill="#89043d") +
    theme_classic() + theme(legend.position="none") + scale_x_discrete(limits=seq(nreps),labels=seq(nreps))
  grid.draw(g)
  dev.off()
}
