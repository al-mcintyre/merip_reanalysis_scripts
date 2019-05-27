library(CLIPanalyze)
library(GenomeInfoDb)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
#library(BSgenome.Mmusculus.UCSC.mm10)
library(BSGenome.Mmusculus.gencode.M19)
#library(GenomicFeatures) #source of select function that gets overwitten by dplyr
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

setwd('/Volumes/Markarian_501/merip_reanalysis/fig2/')
source('plot_fig2_helper.R')

annotation.order <- c("utr3","utr5","exon","intron","utr3*","utr5*")
control.summary <- read.csv('control_summary.txt',sep=' ')
control.summary$studycont <- paste0(control.summary$study,"_",control.summary$control_type)
pos.controls <- control.summary[which(control.summary$control_type == "positive"),]$condition

ncontrols <- dim(control.summary)[1]
tools <- c("DESeq2_GLM","edgeR_GLM","QNB","MeTDiff")
more.tools <- c(tools,c("intersect","union"))
output.columns <- c("tool","experiment","control_type","n_peaks","n_peaks_sig","percent_under0.05","n_reps")
tool.results <- as.data.frame(matrix(ncol=length(output.columns), nrow=ncontrols*length(more.tools)))
colnames(tool.results) <- output.columns
tool.results$tool <- sort(rep(more.tools,ncontrols))
tool.results$experiment <- rep(control.summary$study,length(more.tools))
tool.results$control_type <- rep(control.summary$control_type,length(more.tools))
tool.results$Condition <- rep(control.summary$condition,length(more.tools))

tool.colours <- data.frame(colour=c("#66D7D1","#1b9aaa","#db7b78","#323253"))
row.names(tool.colours) <- tools
#all.counts <- list() 

control.colours <- c(get_palette(c("#ff9b3c","#ffd07a"),length(which(control.summary$control_type == "negative"))),
                     get_palette(c("#7a003b","#dba2ba"),length(which(control.summary$control_type == "positive"))))
c.colours <- data.frame(colour=control.colours)
row.names(c.colours) <- control.summary$studycont

for (expcont in control.summary$studycont[1]){
  expcont <- 'engel_neuron_2018_negative'
  #exp <- control.summary[which(control.summary$studycont == expcont),'study']
  ctype <- control.summary[which(control.summary$studycont == expcont),'control_type']
  col.data <- read.csv(paste0('samplesheet_',ctype,'_',exp,'.csv'))
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
  
  peaks <- import.peaks(unique(col.data$Peaks))
  peak.counts <- get.counts(col.data,peaks,bsg=bsgenome)
  all.counts[[paste0(exp,"_",ctype)]] <- peak.counts
  #all.counts[[expcont]]$peaks <- peaks$peaks
  
  #writeLines("getting counts for full genes")
  gene.counts <- get.gene.counts(col.data[which(col.data$IP == 'input'),],gtf)
  all.counts[[paste0(exp,"_",ctype)]]$genes <- gene.counts
  remove(peaks)
}

plot_overdispersion <- function(means,vars,label){
  png(paste0('sfig2a_',label,"_mean_variance.png"),width=100,height=100,units='mm',res=500)
  plot(vars[which(means != 0 & vars != 0)]~means[which(means != 0 & vars != 0)],
       log="xy",xlab="mean",ylab="variance",cex=0.8,
       pch=16,col = rgb(red = 0, green = 0, blue = 0, alpha = 0.3), mar=c(0,0,0,0))
  reg <- lm(log10(vars[which(means != 0 & vars != 0)])~log10(means[which(means != 0 & vars != 0)]))
  abline(reg,col="#d96c06")
  text(200,5500000,"regression",col="#d96c06")
  abline(0,1,col="#3a5da3")
  text(2000,400,"x = y",col="#3a5da3")
  dev.off()
}

for (exp in control.summary[which(control.summary$control_type == 'negative'),'study']){
  exp <- control.summary[which(control.summary$control_type == 'negative'),'study'][2]
  ctype <- 'negative'
  expcont <- paste0(exp,"_",ctype)
  col.data <- read.csv(paste0('samplesheet_',ctype,'_',exp,'.csv'))
  counts <- counts(all.counts[[paste0(exp,"_",ctype)]]$peak.counts)
  
  #plot mock mean vs. variance for Figure S2A
  meanVar <- data.frame("mean"=rowMeans(counts[,which(col.data$IP == "IP")]),"variance"=rowVar(counts[,which(col.data$IP == "IP")]))
  plot_overdispersion(meanVar$mean,meanVar$variance,exp)
  
  if (exp == "neg_control_mocks"){
    col.data$Annotation <- paste0(col.data$Factor,', ',col.data$Condition)
    pca.col <- c("#e77728","#7dbbc3","#564256",'#947cc1')} else{
      col.data$Annotation <- as.factor(col.data$Replicate)
      pca.col <- c('#987284','#75b9be','#8bb174','#f9b5ac','#ee7674','#426b69','#a74482')}
  
  #plot mock sample QC for Figure S2b
  pdf(paste0('sfig2b_',exp,'_peak_pca.pdf'),height=3.5,width=4.5)
  pca <- autoplot(prcomp(t(counts)), data = col.data, colour = 'Annotation', shape='IP')
  pca + scale_fill_manual(values = pca.col) + scale_color_manual(values = pca.col) + theme(panel.background = element_blank())
  dev.off()
  
  #plot mock sample for Figure S1b  
  mean.inputs4peaks <- rowMedians(counts[,which(col.data$IP =="input")])
  pdf(paste0('sfig1_',exp,'_peaks_by_inputs.pdf'),height=4,width=4)
  qplot(mean.inputs4peaks, geom="histogram", binwidth = 0.03) + scale_x_continuous(trans='log10') + xlab("log10(median input read count)") +
    ylab("# peaks") + theme(panel.background = element_blank())
  dev.off()
  # 
  # if (exp == "engel_neuron_2018"){ rep.cols <- c(2,5)} else{rep.cols <- c(2,10)}
  # ip.over.input <- counts[,which(col.data$IP =="IP")[rep.cols]]/counts[,which(col.data$IP =="input")[rep.cols]]
  # ip.over.input <- as.data.frame(ip.over.input[complete.cases(ip.over.input),])
  # colnames(ip.over.input) <- c("ReplicateA","ReplicateB")
  # tiff(paste0('sfig1_',exp,'_ip_over_input_repsAB.tiff'),units="in",height=3.5,width=4.5,res=500)
  # ggscatter(ip.over.input,x="ReplicateA",y="ReplicateB", add="reg.line", add.params = list(color = "black"),
  #           conf.int = FALSE, cor.coef = TRUE, cor.method = "pearson",
  #           xlab="Replicate A, IP/input (log2)", ylab="Replicate B, IP/input (log2)",color="#89043d",alpha=0.4) +
  #           scale_x_continuous(trans='log2') + scale_y_continuous(trans='log2') 
  # dev.off()
  
  ##Figure 2A, comparison of negative binomial and Poisson model fits using log likelihoods
  
  #compare common dispersion estimates between IP and input (biological coeff of variation)^2
  writeLines(paste('edgeR common dispersion estimate for combined input + IP samples =',estimateCommonDisp(counts)))
  writeLines(paste('edgeR common dispersion estimate for IP samples =',estimateCommonDisp(counts[,which(col.data$IP == "IP")])))
  writeLines(paste('edgeR common dispersion estimate for input samples =',estimateCommonDisp(counts[,which(col.data$IP == "input")])))

  #for Inputs + IPs,
  #for MLE NB fit and MLE Poisson fit generate 500 simulated datasets based on means and dispersions (also tried for DESeq2 NB fit, edgeR NB fit)
  #and compare log likelihoods to original samples

  for (mod in c("Negative Binomial","Poisson")){
    data <- 'ip_input'
    methods <- c('MLE')
    if (data == 'ip_input'){
      subcounts <- counts
      submeta <- col.data
    } else{
      subcounts <- counts[,which(col.data$IP == data)]
      submeta <- col.data[which(col.data$IP == data),]
    }

    mle.md <- apply(unname(subcounts),1,fit.mle,model=mod)
    if (mod == 'Negative Binomial'){
      mle.md <- do.call(rbind, mle.md)
      keep.peaks <- which(complete.cases(mle.md))   #remove peaks where MLE unable to determine fit (18 peaks for mocks)
    } else{ mle.md <- as.data.frame(mle.md)}

    mle.md <- mle.md[keep.peaks,] #keep same peaks for Poisson as Negative Binomial
    subcounts <- subcounts[keep.peaks,]
    mds <- data.frame("mle"=mle.md)

    nexps = dim(subcounts)[2]
    npeaks <- dim(mle.md)[1]
    iterations <- 500
    sim.ll.mat <- matrix(ncol=length(methods), nrow=iterations)
    sam.lls <- list()
    for (ind in 1:length(methods)){
      sim.lls <- list()
      if (mod == "Negative Binomial"){ md <- mds[,c(2*ind-1,2*ind)] } else{ md <- mds}

      for (i in c(1:iterations)){
        simcounts <- t(apply(md,1,sim.data,model=mod))
        simcounts.md <- cbind(simcounts,md)
        nb.ll.sim <- apply(simcounts.md,1,calc.ll,model=mod)
        nb.sum.ll.sim <- mean(nb.ll.sim)
        sim.lls[[i]] <- nb.sum.ll.sim
      }
      sim.ll.mat[,ind] <- unlist(sim.lls)
      counts.md <- cbind(subcounts,md)
      nb.ll.sample <- apply(counts.md,1,calc.ll,model=mod)
      sam.lls[ind] <- mean(nb.ll.sample)
    }

    sam.lls <- unlist(sam.lls)
    writeLines(paste('For',mod,'sample falls into',length(sim.ll.mat[,1][which(sim.ll.mat[,1] < sam.lls[1])])/(iterations/100),'percentile'))
    for (ind in 1:length(methods)){
      method <- methods[ind]
      pdf(paste0('fig2a_ll_',exp,'_',mod,'_fit.pdf'),height=3,width=6)
      par(mar=c(par('mar')[1:3], 0))
      nb.sum.ll.sample <- sam.lls[ind]
      hist(sim.ll.mat[,ind],xlim=c(min(c(sim.ll.mat[,ind],nb.sum.ll.sample)-0.05),max(c(sim.ll.mat[,ind],nb.sum.ll.sample))+0.05),
           xlab="Mean log likelihood\nacross peaks",main="",breaks=20,col='black')
      abline(v=nb.sum.ll.sample,col=as.character(c.colours[expcont,]),lwd = 2)
      legend("topright",legend = c("sample","simulated data"),col = c(as.character(c.colours[expcont,]),"black"),lwd = 2, inset=c(0,0),
             xpd=TRUE,bty='n')
      dev.off()
    }
  }
}


source('run_tools_test.R')

##compare % peaks with p value < 0.05 in all mocks, 3 reps of mocks, DAA data, FTO data, and WTAP data -- add GLM with covariate for non-peak gene coverage

for (expcont in control.summary$studycont[c(1:10)]){
  
  ctype <- control.summary[which(control.summary$studycont == expcont),'control_type']
  exp <- control.summary[which(control.summary$studycont == expcont),'study']
  col.data <- read.csv(paste0('samplesheet_',ctype,'_',exp,'.csv'))
  counts <- counts(all.counts[[paste0(exp,"_",ctype)]]$peak.counts)
  reference_condition = col.data$Condition[1]
  p.results <- list()
  
  for (tool in tools){
    p.adj.l2fc <- run.tool(tool,counts,col.data,reference_condition,full.results = TRUE)
    all.counts[[expcont]][[tool]] <- p.adj.l2fc
    p <- p.adj.l2fc$p
    p.results[[tool]] <- p
    df.inds <- which(tool.results$tool == tool & tool.results$experiment == exp & tool.results$control_type == ctype)
    tool.results[df.inds,"n_peaks"] <- length(p)
    tool.results[df.inds,"n_peaks_sig"] <- length(which(p < 0.05))
    tool.results[df.inds,"percent_under0.05"] <- length(which(p < 0.05))*100/length(p)
    pdf(paste0("pvaluedists/",tool,"_",exp,"_",ctype,"_pvalues.pdf"),width=4,height=2.5)
    hist(p,col=as.character(tool.colours[tool,"colour"]),xlab=paste(tool,"p values"),main="")
    dev.off()
    writeLines(paste(exp,tool,length(which(p.results[[tool]]<0.05) )*100/length(p)))
  }
  p.results[["intersect"]] <- intersect(intersect(which(p.results[["DESeq2_GLM"]]<0.05),which(p.results[["edgeR_GLM"]]<0.05)),which(p.results[["QNB"]]<0.05))
  p.results[["union"]] <- union(union(which(p.results[["DESeq2_GLM"]]<0.05),which(p.results[["edgeR_GLM"]]<0.05)),which(p.results[["QNB"]]<0.05))
  for (tool in c("intersect","union")){
    df.inds <- which(tool.results$tool == tool & tool.results$experiment == exp & tool.results$control_type == ctype)
    writeLines(paste(exp,tool,length(p.results[[tool]])*100/length(p)))
    tool.results[df.inds,"n_peaks"] <- length(p)
    tool.results[df.inds,"n_peaks_sig"] <- length(p.results[[tool]])
    tool.results[df.inds,"percent_under0.05"] <- length(p.results[[tool]])*100/length(p)
  }
  
}

gvp.cols <- c('peak.L2FC','gene.L2FC','padj','tool','experiment')
gene.vs.peak.by.tool <- as.data.frame(matrix(ncol=length(gvp.cols), nrow=0))
pvg.cols <- c('peak.count','tool','gene.sig')
peaks.vs.gene.sig <- as.data.frame(matrix(ncol=length(pvg.cols)),nrow=0)
colnames(gene.vs.peak.by.tool) <- gvp.cols
colnames(peaks.vs.gene.sig) <- pvg.cols

for (expcont in control.summary$studycont[c(3:10)]){
  
  ctype <- 'positive'
  exp <- control.summary[which(control.summary$studycont == expcont),'study']
  col.data <- read.csv(paste0('samplesheet_',ctype,'_',exp,'.csv'))
  ip.counts <- counts(all.counts[[expcont]]$peak.counts)[,which(col.data$IP == 'IP')]
  col.data <- col.data[which(col.data$IP == 'input'),]
  peaks <- all.counts[[expcont]]$peaks
  gene.counts <- all.counts[[expcont]]$genes  
  cont <- as.character(control.summary[which(control.summary$studycont == expcont),'control'])
  cond <- as.character(control.summary[which(control.summary$studycont == expcont),'stimulus'])
  condition.label <- as.character(control.summary[which(control.summary$studycont == expcont),'condition'])
  peak2gene.inds <- sapply(peaks$geneid, function(p) which(rownames(gene.counts) == p)[1])
  
  gene.de <- run.deseq.gene(cond,cont,gene.counts,col.data,"gene")
  ip.de <- run.deseq.gene(cond,cont,ip.counts,col.data,"ip")
  all.counts[[expcont]]$peaks$gene_l2fc <- gene.de$gene_l2fc[peak2gene.inds]
  all.counts[[expcont]]$peaks$gene_p <- gene.de$gene_p[peak2gene.inds]
  all.counts[[expcont]]$peaks$gene_padj <- gene.de$gene_padj[peak2gene.inds]
  all.counts[[expcont]]$peaks$ip_l2fc <- ip.de$ip_l2fc
  all.counts[[expcont]]$peaks$d.l2fc <- all.counts[[expcont]]$peaks$ip_l2fc - all.counts[[expcont]]$peaks$gene_l2fc

  peaks <- all.counts[[expcont]]$peaks
  for (tool in tools){
    sig.peaks <- which(all.counts[[expcont]][[tool]]$padj < 0.05)
    if (length(sig.peaks) > 0){
      sig.peak.df <- as.data.frame(cbind(peaks[sig.peaks,]$ip_l2fc,peaks[sig.peaks,]$gene_l2fc,peaks[sig.peaks,]$gene_padj,tool,condition.label))
      colnames(sig.peak.df) <- gvp.cols
      gene.vs.peak.by.tool <- rbind(gene.vs.peak.by.tool,sig.peak.df)
    }
  }
  
}

summary <- doBy::summaryBy(percent_under0.05 ~ tool + control_type + Condition, tool.results[complete.cases(tool.results[,1:4]),], FUN = c(mean, sd))
summary$percent_under0.05.sd[which(is.na(summary$percent_under0.05.sd))] <- 0
study.medians <- doBy::summaryBy(percent_under0.05.mean ~ control_type + Condition, data=summary, FUN=list(median))
study.order <- study.medians$Condition[c(2,1,order(study.medians[which(study.medians$Condition %in% pos.controls),]$percent_under0.05.mean.median)+2)]
summary$Condition <- factor(summary$Condition,levels = study.order) #unique(control.summary$condition[which(control.summary$condition != 'Mettl3cKO')]))

main.summary <- summary[which(!summary$tool %in% c('union','intersect')),]
pdf("fig2b_controls_percent_changes.pdf",width = 7, height=4)
ggplot(data=main.summary, aes(x=tool, y=percent_under0.05.mean, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge()) + ylab("% peaks p < 0.05") +
  geom_errorbar(data = main.summary, mapping = aes(ymin = percent_under0.05.mean - percent_under0.05.sd, ymax = percent_under0.05.mean + percent_under0.05.sd),position=position_dodge(.9),width=0) +
  xlab("") + geom_hline(yintercept=5,linetype="dashed",color="#800000") +
  scale_fill_manual(values=control.colours) + theme(panel.grid.major.y = element_line(colour = "lightgrey"), panel.background = element_blank())
dev.off()

comb.summary <- summary[which(summary$tool %in% c('union','intersect')),]
pdf("fig2c_controls_tool_combinations.pdf",width = 5, height=4)
ggplot(data=comb.summary, aes(x=tool, y=percent_under0.05.mean, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge()) + ylab("% peaks") +
  geom_errorbar(data = comb.summary, mapping = aes(ymin = percent_under0.05.mean - percent_under0.05.sd, ymax = percent_under0.05.mean + percent_under0.05.sd),position=position_dodge(.9),width=0) +
  xlab("") + geom_hline(yintercept=5,linetype="dashed",color="#800000") +
  scale_fill_manual(values=control.colours) + theme(panel.grid.major.y = element_line(colour = "lightgrey"), panel.background = element_blank())
dev.off()

gene.vs.peak.by.tool$PmG.L2FC <- as.numeric(as.character(gene.vs.peak.by.tool$peak.L2FC)) - as.numeric(as.character(gene.vs.peak.by.tool$gene.L2FC))
gvp <- gene.vs.peak.by.tool[complete.cases(gene.vs.peak.by.tool),]
gvp$peak.L2FC <- as.numeric(as.character(gvp$peak.L2FC))
gvp$gene.L2FC <- as.numeric(as.character(gvp$gene.L2FC))
gvp$tool <- factor(gvp$tool,levels=rev(tools))
pdf(paste0('fig2_peak_vs_geneL2FC.pdf'),height=4,width=4)
ggscatter(gvp,x="peak.L2FC",y="gene.L2FC", add="reg.line", add.params = list(color = "tool"),
          conf.int = FALSE,xlab="Peak L2FC", ylab="Gene L2FC",color="tool",alpha=0.4,palette=as.character(rev(tool.colours$colour)),
          legend.title = '',legend='top',font.legend = c(8, "plain", "black"))
dev.off()  

for (tool in tools){
  subgvp <- gvp[which(gvp$tool == tool),]
  pearson <- cor.test(subgvp$peak.L2FC,subgvp$gene.L2FC)
  writeLines(paste(tool,"peak vs. gene change Pearson's R =",pearson$estimate,' p =',pearson$p.value))
}

ggplot(pvg,aes(x = tool, y = peak.count,fill = gene.sig)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = percent_format())

gvp$tool <- factor(gvp$tool,levels=tools)
gvp$PmG.L2FC <- as.numeric(as.character(gvp$PmG.L2FC))
pdf(paste0('sfig2_peak_minus_geneL2FC.pdf'),height=4,width=8)
ggplot(data=gvp, aes(x=tool,y=PmG.L2FC,fill=tool)) + geom_violin(alpha=0.2) + geom_boxplot(width=0.1,alpha=0.4) +
  xlab("method") + ylab("Peak - gene L2FC") + scale_fill_manual(values=as.character(tool.colours$colour)) +
  theme(panel.background = element_blank(),axis.text.x = element_text(angle=90, hjust=1)) +
  theme_classic() + theme(legend.position="none") + geom_hline(yintercept=0,linetype="dashed",color="#800000") +
  geom_hline(yintercept=1,linetype="dashed",color="pink") + geom_hline(yintercept=-1,linetype="dashed",color="pink") +
  stat_summary(fun.data = function(x) data.frame(y=-8, label = paste(round(median(x),3))), geom="text")
dev.off()


#Figure to show # significant peaks increased and decreased in positive controls 
#this may not be the best measure anyway - eg. maybe enrichment of METTL16-catalyzed sites is increased with k/d of METTL3? Because more antibody available? 
peaks.updown <- gvp %>% group_by(experiment,tool) %>% 
  summarise(down= sum(PmG.L2FC<0), up= sum(PmG.L2FC >=0))
peaks.updown <- reshape2::melt(peaks.updown,id=c("experiment","tool"))
peaks.updown[which(peaks.updown$variable == 'down'),"value"] <- -log10(peaks.updown[which(peaks.updown$variable == 'down'),"value"])
peaks.updown[which(peaks.updown$variable == 'up'),"value"] <- log10(peaks.updown[which(peaks.updown$variable == 'up'),"value"])
peaks.updown[which(is.infinite(peaks.updown$value)),"value"] <- 0
peaks.updown$experiment <- factor(peaks.updown$experiment,levels=study.order[3:length(study.order)])
peaks.updown <- tidyr::complete(peaks.updown, experiment, tool, variable, fill = list(value = 0))

pdf('sfig2_maybe_peaks_updown.pdf',height=4,width=6)
ggplot(peaks.updown, aes(x=tool, y=value, fill=experiment)) +
  geom_bar(stat="identity", position=position_dodge()) + ylab('log10(# peaks up or down)') +
  scale_fill_manual(values=control.colours[3:length(study.order)]) + 
  theme(panel.grid.major.y = element_line(colour = "lightgrey"), panel.background = element_blank()) +
  geom_hline(yintercept=0,color="lightgrey") 
dev.off()