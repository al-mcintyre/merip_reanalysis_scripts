library(CLIPanalyze)
library(GenomeInfoDb)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
#library(BSgenome.Mmusculus.UCSC.mm10)
library(BSGenome.Mmusculus.gencode.M19)
library(GenomicFeatures) #source of select function that gets overwitten by dplyr
library(stringr)
library(VennDiagram)
#library(MASS) #also has select function, required for fitdistr
library(edgeR)
library(ggplot2)
library(doBy)
#biocLite("exomePeak")
library(QNB)
library(moments)
library(scatterplot3d)
library(ggpubr) #for ggscatter
library(ggfortify) #needed for autoplot

library(BiocParallel)
register(SerialParam())
library(DESeq2)


setwd('/Volumes/Markarian_501/merip_reanalysis/fig3/')
source('../fig2/plot_fig2_helper.R')
source('../fig2/run_tools_test.R')

all.counts <- list() 

annotation.order <- c("utr3","utr5","exon","intron","utr3*","utr5*")
experiment.summary <- read.csv('experiment_summary.txt',sep=' ')

nexps <- dim(experiment.summary)[1]
tools <- c("DESeq2_GLM","edgeR_GLM","QNB")
more.tools <- c("original",tools,"filtered")
output.columns <- c("tool","experiment","n_peaks","n_peaks_sig","percent_under0.05","n_reps")
tool.results <- as.data.frame(matrix(ncol=length(output.columns), nrow=nexps*length(more.tools)))
colnames(tool.results) <- output.columns
tool.results$tool <- sort(rep(more.tools,nexps))
tool.results$experiment <- rep(experiment.summary$study,length(more.tools))
tool.results$n_reps <- rep(experiment.summary$reps,length(more.tools))
tool.results$Condition <- rep(experiment.summary$condition,length(more.tools))

tool.cols <- c("#707070","#66D7D1","#1b9aaa","#db7b78","#721652")
tool.colours <- data.frame(colour=tool.cols)
row.names(tool.colours) <- more.tools

for (exp in experiment.summary$study[1:12]){
  col.data <- read.csv(paste0('samplesheet_',exp,'.csv'))
  genomeTag <- as.character(experiment.summary[which(experiment.summary$study == exp),'species'])
  if (genomeTag == "mm10"){
    namespace <- loadNamespace("BSGenome.Mmusculus.gencode.M19")
    bsgenome <- namespace[["BSGenome.Mmusculus.gencode.M19"]]
    gtf <- "/Volumes/Markarian_501/merip_reanalysis/fig2/gencode.vM19.annotation.gtf"
  }else {
    bsgenome <- load.bsgenome(genomeTag)
    gtf <- "/Volumes/Markarian_501/20180228_flaviviridae/gencode.v25.annotation.gtf"
  }
  annot <- loadAnnot(genomeTag)

  ##Section to comment out once all.counts has been filled in (if want to adjust thresholds)
  #peaks <- import.peaks(unique(col.data$Peaks))
  #peak.counts <- get.counts(col.data,peaks,bsg=bsgenome)
  #all.counts[[exp]] <- peak.counts
  #all.counts[[exp]]$peaks <- peaks$peaks
  #peaks <- all.counts[[exp]]
  
  counts <- counts(all.counts[[exp]]$peak.counts)
  
  #gene.counts <- get.gene.counts(col.data[which(col.data$IP == 'input'),],gtf)
  #all.counts[[exp]]$genes <- gene.counts
  
  #gene.counts <- all.counts[[exp]]$genes
  #input.data <- col.data[which(col.data$IP == 'input'),]
  #ip.data <- col.data[which(col.data$IP == 'IP'),]
  #cond <- as.character(experiment.summary[which(experiment.summary$study == exp),]$stimulus)
  #cont <- as.character(experiment.summary[which(experiment.summary$study == exp),]$control)
  #peak2gene.inds <- sapply(peaks$peaks$geneid, function(p) which(rownames(gene.counts) == p)[1])
  
  #gene.de <- run.deseq.gene(cond,cont,gene.counts,input.data,"gene")
  #ip.counts <- counts(all.counts[[exp]]$peak.counts)[,which(col.data$IP == 'IP')]
  #ip.de <- run.deseq.gene(cond,cont,ip.counts,ip.data,"ip")
  #all.counts[[exp]]$peaks$gene_l2fc <- gene.de$gene_l2fc[peak2gene.inds]
  #all.counts[[exp]]$peaks$gene_p <- gene.de$gene_p[peak2gene.inds]
  #all.counts[[exp]]$peaks$gene_padj <- gene.de$gene_padj[peak2gene.inds]
  #all.counts[[exp]]$peaks$ip_l2fc <- ip.de$ip_l2fc
  #all.counts[[exp]]$peaks$d.l2fc <- all.counts[[exp]]$peaks$ip_l2fc - all.counts[[exp]]$peaks$gene_l2fc
  
  PmG <- all.counts[[exp]]$peaks$d.l2fc
  col.data <- read.csv(paste0('samplesheet_',exp,'.csv'))

  s.peaks <- c()
  
  reference_condition = col.data$Condition[1]
  
  for (tool in tools){
    p.adj.l2fc <- all.counts[[exp]][[tool]]
    #p.adj.l2fc <- run.tool(tool,counts,col.data,reference_condition,full.results = TRUE)
    
    p <- p.adj.l2fc$p
    padj <- p.adj.l2fc$padj
    #all.counts[[exp]][[tool]] <- p.adj.l2fc

    df.inds <- which(tool.results$tool == tool & tool.results$experiment == exp)
    tool.results[df.inds,"n_peaks"] <- length(p)
    tool.results[df.inds,"n_peaks_sig"] <- length(which(padj < 0.05))
    tool.results[df.inds,"percent_under0.05"] <- length(which(p < 0.05))*100/length(p)
    pdf(paste0("pvaluedists/",tool,"_",exp,"_pvalues.pdf"),width=4,height=2.5)
    hist(p,col=as.character(tool.colours[tool,"colour"]),xlab=paste(tool,"p values"),main="")
    dev.off()
    writeLines(paste(exp,tool,length(which(p<0.05) )*100/length(p),length(which(padj<0.05))))
    sig.peaks.this.tool <- which(p.adj.l2fc$padj < 0.05)
    with.min <- which(abs(PmG) >= 1 & rowMin(counts) >= 10 & !(all.counts[[exp]]$peaks$annot %in% c("intron","intergenic")))
    s.peaks <- union(s.peaks,intersect(sig.peaks.this.tool,with.min))
    if (length(sig.peaks.this.tool) > 0){
      write.table(as.data.frame(all.counts[[exp]]$peaks[paste0("peak",sig.peaks.this.tool),])[,c("seqnames","start","end","maingene","width","strand","annot","d.l2fc")],paste0(exp,"_",tool,"_sig_peaks.bed"),quote=FALSE,sep="\t",row.names = FALSE, col.names = FALSE)
      write.table(unique(as.data.frame(all.counts[[exp]]$peaks[paste0("peak",sig.peaks.this.tool),])$maingene),paste0(exp,"_",tool,"_sig_genes.bed"),quote=FALSE,sep="\t",row.names = FALSE, col.names = FALSE)
    }
  }
  
  tool <- "filtered"
  df.inds <- which(tool.results$tool == tool & tool.results$experiment == exp)
  tool.results[df.inds,"n_peaks_sig"] <- length(s.peaks)
  tool <- "original"
  df.inds <- which(tool.results$tool == tool & tool.results$experiment == exp)
  tool.results[df.inds,"n_peaks_sig"] <- experiment.summary$reported_peaks[which(experiment.summary$study == exp)]
  if (length(s.peaks) > 0){
    write.table(as.data.frame(all.counts[[exp]]$peaks[paste0("peak",s.peaks),])[,c("seqnames","start","end","maingene","width","strand","annot","d.l2fc")],paste0(exp,"_sig_peaks.bed"),quote=FALSE,sep="\t",row.names = FALSE, col.names = FALSE)
  }
}


summary <- tool.results[complete.cases(tool.results[,c('tool','experiment','n_peaks_sig','n_reps','Condition')]),]
filtered.summary <- summary[which(summary$tool == "filtered"),]
condition.order <- filtered.summary[order(filtered.summary$n_peaks_sig, -filtered.summary$n_reps),'Condition']
summary$Condition <- factor(summary$Condition,levels = condition.order)
summary$tool <- factor(summary$tool,levels = more.tools)
summary$n_peaks_sig <- summary$n_peaks_sig + 1
pdf("fig3a_number_sig_peaks.pdf",width = 6, height=4)
ggplot(data=summary, aes(x=Condition, y=n_peaks_sig, fill=tool)) +
  geom_bar(stat="identity", position=position_dodge()) + ylab("# peaks changed + 1 (log10 scale)") + 
  facet_grid(~ Condition, space="free_x", scales="free_x", switch="x") + 
  xlab("") + scale_y_continuous(trans='log10') +
  scale_fill_manual(values=tool.cols) + theme(panel.grid.major.y = element_line(colour = "lightgrey"), panel.grid.major.x = element_blank(),
                                              panel.background = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.background = element_blank(),strip.text.x = element_blank(),panel.spacing = unit(0.2, "lines"))
dev.off()


reanalyses <- c("hiv","kshv","hcmv")
condition.counts <- list()
for (condition in reanalyses[2]){
  experiment.summary <- read.csv(paste0(condition,'_experiment_summary.txt'),sep=' ')
  sig.genes <- list()
  sig.peaks <- list()
  for (exp in experiment.summary$study[4]){
    exp <- as.character(exp)
    col.data <- read.csv(paste0('samplesheet_',exp,'.csv'))
    genomeTag <- as.character(experiment.summary[which(experiment.summary$study == exp),'species'])[1]
    if (genomeTag == "mm10"){
      namespace <- loadNamespace("BSGenome.Mmusculus.gencode.M19")
      bsgenome <- namespace[["BSGenome.Mmusculus.gencode.M19"]]
      gtf <- "/Volumes/Markarian_501/merip_reanalysis/fig2/gencode.vM19.annotation.gtf"
    }else {
      bsgenome <- load.bsgenome(genomeTag)
      gtf <- "/Volumes/Markarian_501/20180228_flaviviridae/gencode.v25.annotation.gtf"
    }
    annot <- loadAnnot(genomeTag)
    
    peaks <- import.peaks(unique(col.data$Peaks))
    peak.counts <- get.counts(col.data,peaks,bsg=bsgenome)
    condition.counts[[exp]] <- peak.counts
    #condition.counts[[exp]]$peaks <- peaks$peaks
    gene.counts <- counts(peak.counts$gene.counts.nopeaks)
    
    peaks <- condition.counts[[exp]]
    counts <- counts(condition.counts[[exp]]$peak.counts)
    
    reference_condition = as.character(col.data$Condition[1])

    gene.counts <- get.gene.counts(col.data[which(col.data$IP == 'input'),],gtf)
    condition.counts[[exp]]$genes <- gene.counts
    input.data <- col.data[which(col.data$IP == 'input'),]
    ip.data <- col.data[which(col.data$IP == 'IP'),]
    cond <- as.character(experiment.summary[which(experiment.summary$study == exp),]$stimulus)
    cont <- as.character(experiment.summary[which(experiment.summary$study == exp),]$control)
    peak2gene.inds <- sapply(peaks$peaks$geneid, function(p) which(rownames(gene.counts) == p)[1])
    
    gene.de <- run.deseq.gene(cond,cont,gene.counts,input.data,"gene")
    ip.counts <- counts(condition.counts[[exp]]$peak.counts)[,which(col.data$IP == 'IP')]
    ip.de <- run.deseq.gene(cond,cont,ip.counts,ip.data,"ip")
    condition.counts[[exp]]$peaks$gene_l2fc <- gene.de$gene_l2fc[peak2gene.inds]
    condition.counts[[exp]]$peaks$gene_p <- gene.de$gene_p[peak2gene.inds]
    condition.counts[[exp]]$peaks$gene_padj <- gene.de$gene_padj[peak2gene.inds]
    condition.counts[[exp]]$peaks$ip_l2fc <- ip.de$ip_l2fc
    condition.counts[[exp]]$peaks$d.l2fc <- condition.counts[[exp]]$peaks$ip_l2fc - condition.counts[[exp]]$peaks$gene_l2fc
    PmG <- condition.counts[[exp]]$peaks$d.l2fc
    col.data <- read.csv(paste0('samplesheet_',exp,'.csv'))
    
    genes <- c()
    s.peaks <- c()
    for (tool in tools){
      p.adj.l2fc <- run.tool(tool,counts,col.data,reference_condition,full.results = TRUE)
      p <- p.adj.l2fc$p
      padj <- p.adj.l2fc$padj
      condition.counts[[exp]][[tool]] <- p.adj.l2fc
      writeLines(paste(exp,tool,length(which(p<0.05) )*100/length(p)))
      genes <- union(genes,peaks$peaks$maingene[which(p.adj.l2fc$padj < 0.05 & abs(PmG) >= 1)])
      s.peaks <- union(s.peaks,which(p.adj.l2fc$padj < 0.05 & abs(PmG) >= 1))
    }
    sig.genes[[exp]] <- genes
    sig.peaks[[exp]] <- condition.counts[[exp]]$peaks[paste0("peak",s.peaks),]
    write.table(as.data.frame(sig.peaks[[exp]])[,c("seqnames","start","end","maingene","width","strand","annot","d.l2fc")],paste0(exp,"_sig_peaks.bed"),quote=FALSE,sep="\t",row.names = FALSE, col.names = FALSE)
  }
  #hcmv.intersect <- intersect(sig.peaks[["winkler_hcmv_2018"]],sig.peaks[["rubio_hcmv_2018"]])
  #write.table(as.data.frame(hcmv.intersect)[,c("seqnames","start","end")],paste0("hcmv_sig_peaks.bed"),quote=FALSE,sep="\t",row.names = FALSE, col.names = FALSE)
  kshv.intersect.1 <- GenomicRanges::intersect(sig.peaks[["hesser_kshv_2018"]],sig.peaks[["tan_kshv_2018"]])
  write.table(as.data.frame(kshv.intersect.1)[,c("seqnames","start","end")],paste0("kshv_sig_peaks.1.bed"),quote=FALSE,sep="\t",row.names = FALSE, col.names = FALSE)
  kshv.intersect.2 <- GenomicRanges::intersect(sig.peaks[["hesser_kshv_2018"]],sig.peaks[["tan_kshv_2018a"]])
  write.table(as.data.frame(kshv.intersect.2)[,c("seqnames","start","end")],paste0("kshv_sig_peaks.2.bed"),quote=FALSE,sep="\t",row.names = FALSE, col.names = FALSE)
  kshv.intersect.3 <- GenomicRanges::intersect(sig.peaks[["hesser_kshv_2018"]],sig.peaks[["tan_kshv_2018b"]])
  write.table(as.data.frame(kshv.intersect.3)[,c("seqnames","start","end")],paste0("kshv_sig_peaks.3.bed"),quote=FALSE,sep="\t",row.names = FALSE, col.names = FALSE)
  
}

##Transcript to genome coordinates (for FOXM1)

library(ensembldb)
library(EnsDb.Hsapiens.v86)

ensdb <- ensDbFromGff(gff="~/Downloads/Hsapiens.hg38.1.gff3")
edbx <- filter(EnsDb.Hsapiens.v86, filter = ~ seq_name == "12") 
rng_tx <- IRanges(start = c(3435), width = c(163), names = c("NM_202002"))

rng_gnm <- transcriptToGenome(rng_tx, edbx) 
