library(S4Vectors)
library(dplyr)
library(deq)
library(doBy)
library(moments)
library(scatterplot3d)
library(ggplot2)
library(ggpubr) #for ggscatter
library(ggfortify) #needed for autoplot

source('R/plot_fig2_helper_revised.R')
peakcaller <- 'metdiff'
if (peakcaller == "metdiff"){outdir <- 'deq_output_metdiff_peaks/'} else{outdir <- 'deq_output/'}

experiment.summary <- read.csv('summary_files/fig2_experiment_summary.txt',sep=' ')
experiment.summary$studycont <- paste0(experiment.summary$control_type,"_",experiment.summary$study)
pos.controls <- experiment.summary[which(experiment.summary$control_type == "positive"),]$condition

#set up output data frame
ncontrols <- dim(experiment.summary)[1]
tools <- c("deseq2","edger","qnb","metdiff")
more.tools <- c(tools,c("intersect","union"))
output.columns <- c("tool","experiment","control_type","n_peaks","n_peaks_sig","percent_under0.05","n_reps")
tool.results <- as.data.frame(matrix(ncol=length(output.columns), nrow=ncontrols*length(more.tools)))
colnames(tool.results) <- output.columns
tool.results$tool <- sort(rep(more.tools,ncontrols))
tool.results$experiment <- rep(experiment.summary$study,length(more.tools))
tool.results$control_type <- rep(experiment.summary$control_type,length(more.tools))
tool.results$Condition <- rep(experiment.summary$condition,length(more.tools))

#generate colour palettes
tool.colours <- data.frame(colour=c("#66D7D1","#1b9aaa","#db7b78","#323253"))
row.names(tool.colours) <- tools
control.colours <- c(get_palette(c("#ff9b3c","#ffd07a"),length(which(experiment.summary$control_type == "negative"))),
                     get_palette(c("#7a003b","#dba2ba"),length(which(experiment.summary$control_type == "positive"))))
c.colours <- data.frame(colour=control.colours)
row.names(c.colours) <- experiment.summary$studycont

#run DESeq2, edgeR, MeTDiff, and QNB for control experiments
for (contexp in experiment.summary$studycont){ 
  exp <- experiment.summary[which(experiment.summary$studycont == contexp),"study"]
  coldata <- read.csv(paste0('fig2/samplesheet_',contexp,'.csv'))
  ctype <- as.character(experiment.summary[which(experiment.summary$studycont == contexp),"control_type"])

  gtf <- get.gtf(experiment.summary[which(experiment.summary$studycont == contexp),"species"])
  readlen <- experiment.summary[which(experiment.summary$studycont == contexp),"readlen"]
  fraglen <- experiment.summary[which(experiment.summary$studycont == contexp),"fraglen"] 
  pe <- experiment.summary[which(experiment.summary$studycont == contexp),"pe"]
  outfi <- paste0(outdir,contexp,"_deq_results.txt")
  
  control <- as.character(experiment.summary[which(experiment.summary$studycont == contexp),"control"])
  treatment <- as.character(experiment.summary[which(experiment.summary$studycont == contexp),"stimulus"])
  if (peakcaller == "metdiff"){
    if (ctype == 'positive'){
      peak.files <- list.files(paste0('data/',exp,'/',control,'_',treatment,'_MeTDiff_output/MeTDiff_output/'),pattern='peak.sorted.bed',full.names=TRUE)}else{ 
        peak.files <- list.files(paste0('data/',exp,'/',control,'_MeTDiff_output/MeTDiff_output/'),pattern='peak.sorted.bed',full.names=TRUE) } 
    }else {
      peak.files <- as.vector(unique(coldata$Peaks))}
  
  control <- unique(coldata$Condition)[1]
  treatment <- unique(coldata$Condition)[2]
  input.bams <- as.vector(coldata$bam[which(coldata$Condition == control & coldata$IP == "input")])
  ip.bams <- as.vector(coldata$bam[which(coldata$Condition == control & coldata$IP == "IP")])
  treatment.input.bams <- as.vector(coldata$bam[which(coldata$Condition == treatment & coldata$IP == "input")])
  treatment.ip.bams <- as.vector(coldata$bam[which(coldata$Condition == treatment & coldata$IP == "IP")])
  
  results <- deq(input.bams,ip.bams,treatment.input.bams,treatment.ip.bams,peak.files,gtf,tool='deqm',paired.end=pe,outfi=outfi,readlen=readlen,fraglen=fraglen)
}

#Process negative control experiment data for Figure 2a-b
for (exp in control.summary[which(control.summary$control_type == 'negative'),'study']){
  ctype <- 'negative'
  expcont <- paste0(exp,"_",ctype)
  col.data <- read.csv(paste0('fig2/samplesheet_',ctype,'_',exp,'.csv'))
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
  
  #compare replicates IP enrichment for SFig1 
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

##compare % peaks with p value < 0.05 and gene vs. peak log2 fold change
gvp.cols <- c('diff.L2FC','peak.L2FC','gene.L2FC','padj','tool','experiment')
gene.vs.peak.by.tool <- as.data.frame(matrix(ncol=length(gvp.cols), nrow=0))
pvg.cols <- c('peak.count','tool','gene.sig')
peaks.vs.gene.sig <- as.data.frame(matrix(ncol=length(pvg.cols)),nrow=0)
colnames(gene.vs.peak.by.tool) <- gvp.cols
colnames(peaks.vs.gene.sig) <- pvg.cols

for (contexp in experiment.summary$studycont){
  
  exp <- experiment.summary$study[which(experiment.summary$studycont == contexp)]
  ctype <- experiment.summary$control_type[which(experiment.summary$studycont == contexp)]
  condition.label <- as.character(experiment.summary[which(experiment.summary$studycont == contexp),'condition'])
  results <- read.csv(paste0(outdir,contexp,"_deq_results.txt"),sep="\t")
  p.results <- list()
  
  for (tool in tools){
    p <- results[,paste0(tool,'.p')]
    p.results[[tool]] <- p
    df.inds <- which(tool.results$tool == tool & tool.results$experiment == exp & tool.results$control_type == ctype)
    tool.results[df.inds,"n_peaks"] <- length(p)
    tool.results[df.inds,"n_peaks_sig"] <- length(which(p < 0.05))
    tool.results[df.inds,"percent_under0.05"] <- length(which(p < 0.05))*100/length(p)
    pdf(paste0("pvaluedists/",tool,"_",exp,"_",ctype,"_pvalues.pdf"),width=4,height=2.5)
    hist(p,col=as.character(tool.colours[tool,"colour"]),xlab=paste(tool,"p values"),main="")
    dev.off()
    writeLines(paste(contexp,tool,length(which(p<0.05) )*100/length(p)))
    
    sig.peaks <- results[which(results[,paste0(tool,'.padj')]< 0.05),] 
    if (dim(sig.peaks)[1] > 0){
      sig.peak.df <- as.data.frame(cbind(sig.peaks[,c("diff.l2fc","peak.l2fc","gene.l2fc",paste0(tool,'.padj'))],tool,condition.label))
      colnames(sig.peak.df) <- gvp.cols
      gene.vs.peak.by.tool <- rbind(gene.vs.peak.by.tool,sig.peak.df)
    }
  }
  p.results[["intersect"]] <- intersect(intersect(which(p.results[["deseq2"]]<0.05),which(p.results[["edger"]]<0.05)),which(p.results[["qnb"]]<0.05))
  p.results[["union"]] <- union(union(which(p.results[["deseq2"]]<0.05),which(p.results[["edger"]]<0.05)),which(p.results[["qnb"]]<0.05))
  for (tool in c("intersect","union")){
    df.inds <- which(tool.results$tool == tool & tool.results$experiment == exp & tool.results$control_type == ctype)
    writeLines(paste(exp,tool,length(p.results[[tool]])*100/length(p)))
    tool.results[df.inds,"n_peaks"] <- length(p)
    tool.results[df.inds,"n_peaks_sig"] <- length(p.results[[tool]])
    tool.results[df.inds,"percent_under0.05"] <- length(p.results[[tool]])*100/length(p)
  }
}

summary <- doBy::summaryBy(percent_under0.05 ~ tool + control_type + Condition, tool.results[complete.cases(tool.results[,1:4]),], FUN = c(mean, sd))
summary$percent_under0.05.sd[which(is.na(summary$percent_under0.05.sd))] <- 0
study.medians <- doBy::summaryBy(percent_under0.05.mean ~ control_type + Condition, data=summary, FUN=list(median))
study.order <- study.medians$Condition[c(2,1,order(study.medians[which(study.medians$Condition %in% pos.controls),]$percent_under0.05.mean.median)+2)]
summary$Condition <- factor(summary$Condition,levels = unique(study.order)) #unique(control.summary$condition[which(control.summary$condition != 'Mettl3cKO')]))

main.summary <- summary[which(!summary$tool %in% c('union','intersect')),]
pdf(paste0("fig2b_controls_percent_changes_",peakcaller,".pdf"),width = 7, height=4)
ggplot(data=main.summary, aes(x=tool, y=percent_under0.05.mean, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge()) + ylab("% peaks p < 0.05") +
  geom_errorbar(data = main.summary, mapping = aes(ymin = percent_under0.05.mean - percent_under0.05.sd, ymax = percent_under0.05.mean + percent_under0.05.sd),position=position_dodge(.9),width=0) +
  xlab("") + geom_hline(yintercept=5,linetype="dashed",color="#800000") +
  scale_fill_manual(values=control.colours) + theme(panel.grid.major.y = element_line(colour = "lightgrey"), panel.background = element_blank())
dev.off()

comb.summary <- summary[which(summary$tool %in% c('union','intersect')),]
pdf(paste0("fig2c_controls_tool_combinations_",peakcaller,".pdf"),width = 5, height=4)
ggplot(data=comb.summary, aes(x=tool, y=percent_under0.05.mean, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge()) + ylab("% peaks") +
  geom_errorbar(data = comb.summary, mapping = aes(ymin = percent_under0.05.mean - percent_under0.05.sd, ymax = percent_under0.05.mean + percent_under0.05.sd),position=position_dodge(.9),width=0) +
  xlab("") + geom_hline(yintercept=5,linetype="dashed",color="#800000") +
  scale_fill_manual(values=control.colours) + theme(panel.grid.major.y = element_line(colour = "lightgrey"), panel.background = element_blank())
dev.off()

gvp <- gene.vs.peak.by.tool[complete.cases(gene.vs.peak.by.tool),]
gvp$peak.L2FC <- as.numeric(as.character(gvp$peak.L2FC))
gvp$gene.L2FC <- as.numeric(as.character(gvp$gene.L2FC))
gvp$tool <- factor(gvp$tool,levels=rev(tools))
pdf(paste0('fig2_peak_vs_geneL2FC_',peakcaller,'.pdf'),height=4,width=4)
ggscatter(gvp,x="peak.L2FC",y="gene.L2FC", add="reg.line", add.params = list(color = "tool"),
          conf.int = FALSE,xlab="Peak L2FC", ylab="Gene L2FC",color="tool",alpha=0.4,palette=as.character(rev(tool.colours$colour)),
          legend.title = '',legend='top',font.legend = c(8, "plain", "black"))
dev.off()  

for (tool in tools){
  subgvp <- gvp[which(gvp$tool == tool),]
  pearson <- cor.test(subgvp$peak.L2FC,subgvp$gene.L2FC)
  writeLines(paste(tool,"peak vs. gene change Pearson's R =",pearson$estimate,' p =',pearson$p.value))
}

gvp$tool <- factor(gvp$tool,levels=tools)
gvp$diff.L2FC <- as.numeric(as.character(gvp$diff.L2FC))
pdf(paste0('sfig2_peak_minus_geneL2FC_',peakcaller,'.pdf'),height=4,width=8)
ggplot(data=gvp, aes(x=tool,y=diff.L2FC,fill=tool)) + geom_violin(alpha=0.2) + geom_boxplot(width=0.1,alpha=0.4) +
  xlab("method") + ylab("Peak - gene L2FC") + scale_fill_manual(values=as.character(tool.colours$colour)) +
  theme(panel.background = element_blank(),axis.text.x = element_text(angle=90, hjust=1)) +
  theme_classic() + theme(legend.position="none") + geom_hline(yintercept=0,linetype="dashed",color="#800000") +
  geom_hline(yintercept=1,linetype="dashed",color="pink") + geom_hline(yintercept=-1,linetype="dashed",color="pink") +
  stat_summary(fun.data = function(x) data.frame(y=-8, label = paste(round(median(x),3))), geom="text")
dev.off()


#Figure to show # significant peaks increased and decreased in positive controls 
#this may not be the best measure anyway - eg. maybe enrichment of METTL16-catalyzed sites is increased with k/d of METTL3? Because more antibody available? 
peaks.updown <- gvp %>% group_by(experiment,tool) %>% 
  dplyr::summarise(down= sum(diff.L2FC<0), up= sum(diff.L2FC >=0))
peaks.updown <- reshape2::melt(peaks.updown,id=c("experiment","tool"))
peaks.updown[which(peaks.updown$variable == 'down'),"value"] <- -log10(peaks.updown[which(peaks.updown$variable == 'down'),"value"])
peaks.updown[which(peaks.updown$variable == 'up'),"value"] <- log10(peaks.updown[which(peaks.updown$variable == 'up'),"value"])
peaks.updown[which(is.infinite(peaks.updown$value)),"value"] <- 0
peaks.updown$experiment <- factor(peaks.updown$experiment,levels=study.order[3:length(study.order)])
peaks.updown <- tidyr::complete(peaks.updown, experiment, tool, variable, fill = list(value = 0))

pdf(paste0('sfig2_maybe_peaks_updown_',peakcaller,'.pdf'),height=4,width=6)
ggplot(peaks.updown, aes(x=tool, y=value, fill=experiment)) +
  geom_bar(stat="identity", position=position_dodge()) + ylab('log10(# peaks up or down)') +
  scale_fill_manual(values=control.colours[3:length(study.order)]) + 
  theme(panel.grid.major.y = element_line(colour = "lightgrey"), panel.background = element_blank()) +
  geom_hline(yintercept=0,color="lightgrey") 
dev.off()
