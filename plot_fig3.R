library(S4Vectors)
library(dplyr)
library(deq)
library(doBy)
library(moments)
library(scatterplot3d)
library(ggplot2)
library(ggpubr) #for ggscatter
library(ggfortify) #needed for autoplot

args <- commandArgs(TRUE)
peakcaller <- args[1]
exp.summary <- args[2]
mm10.gtf <- args[3]
hg38.gtf <- args[4]
run.deq <- args[5]

source('plot_fig3_helper.R')
outdir <- paste0('fig3/deq_output_',peakcaller,'/')
dir.create(file.path('fig3',paste0('deq_output_',peakcaller)), showWarnings = FALSE)

experiment.summary <- read.csv(exp.summary,sep=' ')
experiment.summary <- experiment.summary[which(experiment.summary$fig == "3"),]
experiment.summary$studycont <- paste0(experiment.summary$control_type,"_",experiment.summary$study)
pos.controls <- experiment.summary[which(experiment.summary$control_type == "positive"),]$stimulus
set.seed(147)

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
tool.results$Condition <- rep(experiment.summary$stimulus,length(more.tools))

#generate colour palettes
tool.colours <- data.frame(colour=c("#66D7D1","#1b9aaa","#db7b78","#323253"))
row.names(tool.colours) <- tools
control.colours <- c(get_palette(c("#ff9b3c","#ffd07a"),length(which(experiment.summary$control_type == "negative"))),
                     get_palette(c("#7a003b","#dba2ba"),length(which(experiment.summary$control_type == "positive"))))
c.colours <- data.frame(colour=control.colours)
row.names(c.colours) <- experiment.summary$studycont

#run DESeq2, edgeR, MeTDiff, and QNB for control experiments
if (run.deq == "TRUE"){
for (contexp in experiment.summary$studycont){ 
  subexpsum <- experiment.summary[which(experiment.summary$studycont == contexp),]
  exp <- subexpsum$study
  ctype <- as.character(experiment.summary[which(experiment.summary$studycont == contexp),"control_type"])

  gtf <- get.gtf(subexpsum$species,mm10.gtf,hg38.gtf)
  readlen <- subexpsum$readlen 
  fraglen <- subexpsum$fraglen 
  pe <- subexpsum$pe 
  outfi <- paste0(outdir,contexp,"_deq_results.txt")
  
  control <- as.character(subexpsum$control)
  stimulus <- as.character(subexpsum$stimulus)
  
  if (ctype == "positive"){
    input.bams <- list.files(paste0(exp,'/alignments'),pattern=paste0('^',control,'_input_[1-9].star.sorted.bam$'),full.names=TRUE)
    ip.bams <- list.files(paste0(exp,'/alignments'),pattern=paste0('^',control,'_IP_[1-9].star.sorted.bam$'),full.names=TRUE)
    stimulus.input.bams <-list.files(paste0(exp,'/alignments'),pattern=paste0('^',stimulus,'_input_[1-9].star.sorted.bam$'),full.names=TRUE)
    stimulus.ip.bams <- list.files(paste0(exp,'/alignments'),pattern=paste0('^',stimulus,'_IP_[1-9].star.sorted.bam$'),full.names=TRUE)
    if (peakcaller == "metdiff"){
        peak.files <- list.files(paste0(exp,'/',control,'_',stimulus,'_MeTDiff_output/MeTDiff_output'),pattern='peak.sorted.bed',full.names=TRUE)
    }else{
        peak.files <- c(list.files(paste0(exp,'/macs2_results'),pattern=paste0("^",control,'_[1-9].macs2_peaks.narrowPeak.peaks.bed'),full.names=TRUE),
                                                            list.files(paste0(exp,'/macs2_results'),pattern=paste0("^",stimulus,'_[1-9].macs2_peaks.narrowPeak.peaks.bed'),full.names=TRUE))
    }
  }else{ #moving away from sample sheets, but this was the easiest way to deal with the negative controls for now
    coldata <- read.csv(paste0('fig3/samplesheet_',contexp,'.csv'))
    if (peakcaller == "metdiff"){
        peak.files <- list.files(paste0(exp,'/',control,'_MeTDiff_output/MeTDiff_output/'),pattern='peak.sorted.bed',full.names=TRUE)
    }else{
        peak.files <- unique(as.vector(coldata$Peaks))
    }
    control <- as.character(unique(coldata$Condition)[1])
    stimulus <- as.character(unique(coldata$Condition)[2])
    input.bams <- as.vector(coldata$bam[which(coldata$Condition == control & coldata$IP == "input")])
    ip.bams <- as.vector(coldata$bam[which(coldata$Condition == control & coldata$IP == "IP")])
    stimulus.input.bams <- as.vector(coldata$bam[which(coldata$Condition == stimulus & coldata$IP == "input")])
    stimulus.ip.bams <- as.vector(coldata$bam[which(coldata$Condition == stimulus & coldata$IP == "IP")])
  }
  #write.table(input.bams)
  #write.table(ip.bams)
  #write.table(stimulus.input.bams)
  #write.table(stimulus.ip.bams)
  #write.table(peak.files)
  results <- deq(input.bams,ip.bams,stimulus.input.bams,stimulus.ip.bams,peak.files,gtf,tool='deqm',paired.end=pe,outfi=outfi,readlen=readlen,fraglen=fraglen)
}
}

#Process negative control experiment data for Figure 3a-b
for (contexp in experiment.summary[which(experiment.summary$control_type == 'negative'),'studycont']){
  ctype <- 'negative'
  subexpsum <- experiment.summary[which(experiment.summary$studycont == contexp),]
  exp <- subexpsum$study
  coldata <- read.csv(paste0('fig3/samplesheet_',contexp,'.csv'))
  counts <- read.csv(paste0(outdir,contexp,"_deq_results.counts.txt"),sep="\t") # counts(all.counts[[paste0(exp,"_",ctype)]]$peak.counts)
  write.table(head(counts))
  
  #plot mock mean vs. variance for Figure S3A
  meanVar <- data.frame("mean"=rowMeans(counts[,which(coldata$IP == "IP")]),"variance"=rowVar(counts[,which(coldata$IP == "IP")]))
  plot_overdispersion(meanVar$mean,meanVar$variance,exp)
  
  if (exp == "neg_control_mocks"){
    coldata$Annotation <- paste0(coldata$Factor,', ',coldata$Condition)
    pca.col <- c("#e77728","#7dbbc3","#564256",'#947cc1')} else{
      coldata$Annotation <- as.factor(coldata$Replicate)
      pca.col <- c('#987284','#75b9be','#8bb174','#f9b5ac','#ee7674','#426b69','#a74482')}
  
  #plot mock sample QC for Figure S3b
  pdf(paste0('fig3/sfig3b_',exp,'_peak_pca.pdf'),height=3.5,width=4.5)
  pca <- autoplot(prcomp(t(counts)), data = coldata, colour = 'Annotation', shape='IP')
  pca + scale_fill_manual(values = pca.col) + scale_color_manual(values = pca.col) + theme(panel.background = element_blank())
  print(pca)
  dev.off()
  
  #plot replicates for Figure S1d  
  #write.table(coldata)
  #write.table(head(counts))
  median.inputs4peaks <- rowMedians(counts[,which(coldata$IP =="input")])
  pdf(paste0('fig1/sfig1d_',exp,'_peaks_by_inputs.pdf'),height=4,width=4)
  g <- qplot(median.inputs4peaks, geom="histogram", binwidth = 0.03) + scale_x_continuous(trans='log10') + xlab("log10(median input read count)") +
    ylab("# peaks") + theme(panel.background = element_blank())
  print(g)
  dev.off()
  
  #compare replicates IP enrichment for SFig1f
  if (exp == "engel_neuron_2018"){ rep.cols <- c(2,5)} else{rep.cols <- c(2,10)}
  ip.over.input <- counts[,which(coldata$IP =="IP")[rep.cols]]/counts[,which(coldata$IP =="input")[rep.cols]]
  ip.over.input <- as.data.frame(ip.over.input[complete.cases(ip.over.input),])
  colnames(ip.over.input) <- c("ReplicateA","ReplicateB")
  pdf(paste0('fig1/sfig1f_',exp,'_ip_over_input_repsAB.pdf'),height=3.5,width=4.5) #tiff'),units="in",height=3.5,width=4.5,res=500)
  g <- ggscatter(ip.over.input,x="ReplicateA",y="ReplicateB", add="reg.line", add.params = list(color = "black"),
             conf.int = FALSE, cor.coef = TRUE, cor.method = "pearson",
             xlab="Replicate A, IP/input (log2)", ylab="Replicate B, IP/input (log2)",color="#89043d",alpha=0.4) +
             scale_x_continuous(trans='log2') + scale_y_continuous(trans='log2') 
  print(g)
  dev.off()
  
  ##Figure 3A, comparison of negative binomial and Poisson model fits using log likelihoods
  
  #compare common dispersion estimates between IP and input (biological coeff of variation)^2
  #writeLines(paste('edgeR common dispersion estimate for combined input + IP samples =',estimateCommonDisp(counts)))
  #writeLines(paste('edgeR common dispersion estimate for IP samples =',estimateCommonDisp(counts[,which(coldata$IP == "IP")])))
  #writeLines(paste('edgeR common dispersion estimate for input samples =',estimateCommonDisp(counts[,which(coldata$IP == "input")])))
  
  #for Inputs + IPs,
  #for MLE NB fit and MLE Poisson fit generate 500 simulated datasets based on means and dispersions (also tried for DESeq2 NB fit, edgeR NB fit)
  #and compare log likelihoods to original samples
  
  for (mod in c("Negative Binomial","Poisson")){
    data <- 'ip_input'
    methods <- c('MLE')
    if (data == 'ip_input'){
      subcounts <- counts
    } else{
      subcounts <- counts[,which(coldata$IP == data)]
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
      pdf(paste0('fig3/fig3a_ll_',exp,'_',mod,'_fit.pdf'),height=3,width=6)
      par(mar=c(par('mar')[1:3], 0))
      nb.sum.ll.sample <- sam.lls[ind]
      hist(sim.ll.mat[,ind],xlim=c(min(c(sim.ll.mat[,ind],nb.sum.ll.sample)-0.05),max(c(sim.ll.mat[,ind],nb.sum.ll.sample))+0.05),
           xlab="Mean log likelihood\nacross peaks",main="",breaks=20,col='black')
      abline(v=nb.sum.ll.sample,col=as.character(c.colours[contexp,]),lwd = 2)
      legend("topright",legend = c("sample","simulated data"),col = c(as.character(c.colours[contexp,]),"black"),lwd = 2, inset=c(0,0),
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
diffs <- as.data.frame(matrix(ncol=2, nrow=0))

for (contexp in experiment.summary$studycont){
  
  exp <- experiment.summary$study[which(experiment.summary$studycont == contexp)]
  ctype <- experiment.summary$control_type[which(experiment.summary$studycont == contexp)]
  condition.label <- as.character(experiment.summary[which(experiment.summary$studycont == contexp),'stimulus'])
  results <- read.csv(paste0(outdir,contexp,"_deq_results.txt"),sep="\t")
  p.results <- list()
  diffs <- rbind(diffs,cbind(condition.label,results$diff.l2fc))
  
  for (tool in tools){
    p <- results[,paste0(tool,'.p')]
    p.results[[tool]] <- p
    df.inds <- which(tool.results$tool == tool & tool.results$experiment == exp & tool.results$control_type == ctype)
    tool.results[df.inds,"n_peaks"] <- length(p)
    tool.results[df.inds,"n_peaks_sig"] <- length(which(p < 0.05))
    tool.results[df.inds,"percent_under0.05"] <- length(which(p < 0.05))*100/length(p)
    pdf(paste0("fig3/pvaluedists/",tool,"_",exp,"_",ctype,"_pvalues.pdf"),width=4,height=2.5)
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
#write.table(study.medians)
study.order <- study.medians$Condition[c(1,2,order(study.medians[which(study.medians$Condition %in% pos.controls),]$percent_under0.05.mean.median)+2)]
write.table(study.order)
summary$Condition <- factor(summary$Condition,levels = unique(study.order)) #unique(control.summary$condition[which(control.summary$condition != 'Mettl3cKO')]))

main.summary <- summary[which(!summary$tool %in% c('union','intersect')),]
pdf(paste0("fig3/fig3b_controls_percent_changes_",peakcaller,".pdf"),width = 7, height=4)
ggplot(data=main.summary, aes(x=tool, y=percent_under0.05.mean, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge()) + ylab("% peaks p < 0.05") +
  geom_errorbar(data = main.summary, mapping = aes(ymin = percent_under0.05.mean - percent_under0.05.sd, ymax = percent_under0.05.mean + percent_under0.05.sd),position=position_dodge(.9),width=0) +
  xlab("") + geom_hline(yintercept=5,linetype="dashed",color="#800000") +
  scale_fill_manual(values=control.colours) + theme(panel.grid.major.y = element_line(colour = "lightgrey"), panel.background = element_blank())
dev.off()

writeLines('c')

comb.summary <- summary[which(summary$tool %in% c('union','intersect')),]
pdf(paste0("fig3/fig3e_controls_tool_combinations_",peakcaller,".pdf"),width = 5, height=4)
ggplot(data=comb.summary, aes(x=tool, y=percent_under0.05.mean, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge()) + ylab("% peaks") +
  geom_errorbar(data = comb.summary, mapping = aes(ymin = percent_under0.05.mean - percent_under0.05.sd, ymax = percent_under0.05.mean + percent_under0.05.sd),position=position_dodge(.9),width=0) +
  xlab("") + geom_hline(yintercept=5,linetype="dashed",color="#800000") +
  scale_fill_manual(values=control.colours) + theme(panel.grid.major.y = element_line(colour = "lightgrey"), panel.background = element_blank())
dev.off()

writeLines('d')

gvp <- gene.vs.peak.by.tool[complete.cases(gene.vs.peak.by.tool),]
gvp$peak.L2FC <- as.numeric(as.character(gvp$peak.L2FC))
gvp$gene.L2FC <- as.numeric(as.character(gvp$gene.L2FC))
gvp$tool <- factor(gvp$tool,levels=rev(tools))
pdf(paste0('fig3/fig3c_peak_vs_geneL2FC_',peakcaller,'.pdf'),height=4,width=4)
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
pdf(paste0('fig3/sfig3e_peak_minus_geneL2FC_',peakcaller,'.pdf'),height=4,width=8)
ggplot(data=gvp, aes(x=tool,y=diff.L2FC,fill=tool)) + geom_violin(alpha=0.2) + geom_boxplot(width=0.1,alpha=0.4) +
  xlab("method") + ylab("Peak - gene L2FC") + scale_fill_manual(values=as.character(tool.colours$colour)) +
  theme(panel.background = element_blank(),axis.text.x = element_text(angle=90, hjust=1)) +
  theme_classic() + theme(legend.position="none") + geom_hline(yintercept=0,linetype="dashed",color="#800000") +
  geom_hline(yintercept=1,linetype="dashed",color="pink") + geom_hline(yintercept=-1,linetype="dashed",color="pink") +
  stat_summary(fun.data = function(x) data.frame(y=-8, label = paste(round(median(x),3))), geom="text")
dev.off()


#Figure to show # significant peaks increased and decreased in positive controls 
#this may not be the best measure - maybe enrichment of METTL16-catalyzed sites is increased with k/d of METTL3? Because more antibody available? 
peaks.updown <- gvp %>% group_by(experiment,tool) %>% 
  dplyr::summarise(down= sum(diff.L2FC<0), up= sum(diff.L2FC >=0))
peaks.updown <- reshape2::melt(peaks.updown,id=c("experiment","tool"))
peaks.updown[which(peaks.updown$variable == 'down'),"value"] <- -log10(peaks.updown[which(peaks.updown$variable == 'down'),"value"])
peaks.updown[which(peaks.updown$variable == 'up'),"value"] <- log10(peaks.updown[which(peaks.updown$variable == 'up'),"value"])
peaks.updown[which(is.infinite(peaks.updown$value)),"value"] <- 0
peaks.updown$experiment <- factor(peaks.updown$experiment,levels=study.order[3:length(study.order)])
peaks.updown <- tidyr::complete(peaks.updown, experiment, tool, variable, fill = list(value = 0))

pdf(paste0('fig3/sfig3_response_peaks_updown_',peakcaller,'.pdf'),height=4,width=6)
ggplot(peaks.updown, aes(x=tool, y=value, fill=experiment)) +
  geom_bar(stat="identity", position=position_dodge()) + ylab('log10(# peaks up or down)') +
  scale_fill_manual(values=control.colours[3:length(study.order)]) + 
  theme(panel.grid.major.y = element_line(colour = "lightgrey"), panel.background = element_blank()) +
  geom_hline(yintercept=0,color="lightgrey") 
dev.off()

colnames(diffs) <- c('condition','diff.l2fc')
diffs$condition <- factor(diffs$condition,levels=study.order)
diffs$diff.l2fc <- as.numeric(as.character(diffs$diff.l2fc))
diffs <- diffs %>% group_by(condition) %>%
  mutate(skew = skewness(diff.l2fc,na.rm=TRUE))
diffs$skew <- factor(as.character(round(diffs$skew,2)))
pdf(paste0('fig3/sfig3c_peak_diff_distributions_',peakcaller,'.pdf'),height=4,width=6)
ggplot(diffs, aes(x=diff.l2fc, fill = condition)) + geom_density(alpha = 1, position="identity")+ 
  scale_fill_manual(values=control.colours)+ facet_wrap( ~condition, ncol=5) + xlim(-2,2) +
  geom_vline(xintercept = 0,colour="darkgrey")
  #geom_text(data=diffs,x=2,y=2, aes(label=skew))
dev.off()

pdf(paste0('fig3/sfig3_peak_diff_distributions_overlap_',peakcaller,'.pdf'),height=4,width=6)
ggplot(diffs, aes(x=diff.l2fc, fill = condition)) + geom_density(alpha = 0.8, position="identity")+ 
  scale_fill_manual(values=control.colours) + xlim(-2,2)
dev.off()
