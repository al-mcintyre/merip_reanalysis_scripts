library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(reshape2)
library(dplyr)

source('R/plot_sfig1a.R')
#change to wherever your mouse and human gtfs are kept
mm10.gtf <- ""
hg38.gtf <- ""

#number of reads vs number of peaks
read.label <- "input"
reads.v.peaks <- read.csv('fig1/reads_vs_peaks.txt',sep=" ")
#reads.v.peaks <- reads.v.peaks[complete.cases(reads.v.peaks),]
reads.v.peaks$label <- gsub('_',' ',reads.v.peaks$label)

pdf(paste0('sfig1d_',read.label,'_reads_vs_peaks.pdf'),height=6,width=7)
ggscatter(reads.v.peaks,x="reads",y="peaks", add="reg.line",add.params = list(color = "black"),
          conf.int = FALSE, cor.coef = TRUE, cor.method = "pearson",
          xlab=paste("#",read.label,"reads"), ylab="# peaks",color="label",shape="label",alpha=1,
          legend.title = '',legend='top',font.legend = c(7, "plain", "black")) + scale_shape_manual(values=rep(c(8,seq(15,18)),8))
dev.off()  

#plot DRACN enrichment
annotation.order <- c("utr3","utr5","exon","intron","utr3*","utr5*")
control.summary <- read.csv('summary_files/fig2_experiment_summary.txt',sep=' ')
control.summary <- control.summary[which(control.summary$control_type == "negative"),]
control.summary$studycont <- paste0(control.summary$study,"_",control.summary$control_type)
output.columns <- c("experiment","replicates","percent.dracn","num.reads")

for (expcont in control.summary$studycont){
  exp <- control.summary[which(control.summary$studycont == expcont),'study']
  ctype <- control.summary[which(control.summary$studycont == expcont),'control_type']
  col.data <- read.csv(paste0('fig2/samplesheet_',ctype,'_',exp,'.csv'))
  if (exp == "neg_control_mocks"){col.data$Factor <- 'mock'}
  genomeTag <- as.character(control.summary[which(control.summary$study == exp & control.summary$control_type == ctype),'species'])
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
  
  label <- unique(col.data$Factor)
  dracn.summary <- as.data.frame(matrix(ncol=length(output.columns), nrow=0))
  colnames(dracn.summary) <- output.columns
  drac.per.len <- as.data.frame(matrix(ncol=3,nrow=0))
  nreps <- max(c(max(col.data$Replicate),length(col.data$Replicate)/2))
  
  for (replicate in seq(nreps)){
    macs2.file <- paste0("data/",exp,"/macs2_results/",label,"_replicate_intersect_",replicate,".bed")
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
  
  writeLines(paste('% DRACN with 1 replicate =',dracn.summary[1,3]))
  writeLines(paste('% DRACN with 2 replicates =',dracn.summary[2,3]))
  writeLines(paste('% DRACN with',nreps,'replicates =',dracn.summary[nreps,3]))          
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
  
  pdf(paste0('sfig1c_',exp,'_macs2_motif_enrichment.pdf'),height=3.2,width=4.3)
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
