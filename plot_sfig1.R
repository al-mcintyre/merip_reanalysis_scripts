library(CLIPanalyze)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(reshape2)

setwd('/Volumes/Markarian_501/merip_reanalysis/fig1b/')
source('../fig2/plot_fig2_helper.R')
source('../fig2/run_tools_test.R')

#number of reads vs number of peaks
read.label <- "input"
reads.v.peaks <- read.csv('reads_vs_peaks.txt',sep=" ")
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
control.summary <- read.csv('../fig2/control_summary.txt',sep=' ')
control.summary <- control.summary[which(control.summary$control_type == "negative"),]
control.summary$studycont <- paste0(control.summary$study,"_",control.summary$control_type)
output.columns <- c("experiment","replicates","percent.dracn","num.reads")

for (expcont in control.summary$studycont[1]){
  exp <- control.summary[which(control.summary$studycont == expcont),'study']
  ctype <- control.summary[which(control.summary$studycont == expcont),'control_type']
  col.data <- read.csv(paste0('../fig2/samplesheet_',ctype,'_',exp,'.csv'))
  if (exp == "neg_control_mocks"){col.data$Factor <- 'mock'}
  genomeTag <- as.character(control.summary[which(control.summary$study == exp & control.summary$control_type == ctype),'species'])
  #annot.name.id <- setNames(as.list(annot$genenames$gene_id), annot$genenames$gene_name)
  if (genomeTag == "mm10"){
    namespace <- loadNamespace("BSGenome.Mmusculus.gencode.M19")
    bsgenome <- namespace[["BSGenome.Mmusculus.gencode.M19"]]
    gtf <- "/Volumes/Markarian_501/merip_reanalysis/fig2/gencode.vM19.annotation.gtf"
  }else {
    bsgenome <- load.bsgenome(genomeTag)
    gtf <- "/Volumes/Markarian_501/20180228_flaviviridae/gencode.v25.annotation.gtf"
  }
  annot <- loadAnnot(genomeTag)
  label <- unique(col.data$Factor)
  dracn.summary <- as.data.frame(matrix(ncol=length(output.columns), nrow=0))
  colnames(dracn.summary) <- output.columns
  drac.per.len <- as.data.frame(matrix(ncol=3,nrow=0))
  nreps <- max(c(max(col.data$Replicate),length(col.data$Replicate)/2))
  
  for (replicate in seq(nreps)){
    macs2.file <- paste0("macs2_peaks/",label,"_replicate_intersect_",replicate,".bed")
    if(file.size(macs2.file) > 0){
      peaks <- import.peaks(macs2.file)
      peaks <- finish.annotation(col.data,peaks,bsg=bsgenome)
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
  
  pdf(paste0('sfig1b_',exp,'_macs2_motif_enrichment.pdf'),height=3.2,width=4.3)
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

#ggplot(data=dracn.summary2, aes(x=replicates)) + geom_bar(stat="identity",aes(y=value,fill = variable)) +
#  scale_fill_manual(values=c("#b2527b","#89043d"),name = "") + xlab("# Replicates") + ylab(NULL) +
#  theme(panel.background = element_blank(),axis.text.x = element_text(angle=90, hjust=1)) +
#facet_wrap(~ variable, ncol=1,scales = "free_y",strip.position = "left", shrink = FALSE,
#           labeller = as_labeller(c(peaks = "Replicate overlap", percent.dracn= "% DRACN") )) +
#  theme_classic() + theme(legend.position="none") + scale_x_discrete(limits=seq(nreps),labels=seq(nreps))
