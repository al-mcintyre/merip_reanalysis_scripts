library(reshape2)
library(ggplot2)
library(ggpubr)

merip_rt_qpcr_results = as.data.frame(read.csv('merip_rtqpcr_summary_by_gene.txt',header=TRUE,sep="\t"))
peaks.gene.assoc = as.data.frame(peak.table[unique(merip_rt_qpcr_results$Peak),c("gene_names","dengue.PminusG.l2fc","hcv.PminusG.l2fc","zikv.PminusG.l2fc","wnv.PminusG.l2fc")])
colnames(peaks.gene.assoc) <- c("Gene","dengue","hcv","zikv","wnv")
peaks.gene.assoc <- melt(peaks.gene.assoc,id="Gene")
colnames(peaks.gene.assoc) <- c("Gene","Virus","MeRIPseq.change")
merip.rt.qpcr.v.seq <- as.data.frame(merge(peaks.gene.assoc, merip_rt_qpcr_results, by=c("Gene","Virus")))
merip.rt.qpcr.v.seq$Gene <- factor(unlist(merip.rt.qpcr.v.seq$Gene))

pdf("sfig4_merip_seq_vs_rtqpcr.pdf",height = 5.2,width=5)
coul <- get_palette(c("#8e1b64","#db7889","#303b49","#a5c7d1"),length(unique(merip.rt.qpcr.v.seq$Gene)))
ggscatter(merip.rt.qpcr.v.seq,x="Mean",y="MeRIPseq.change", color ="Gene", 
          palette=coul,
          add="reg.line", 
          add.params = list(color = "Gene",palette=coul),
          conf.int = FALSE, cor.coef = TRUE, cor.method = "pearson", 
          xlab="MeRIP-RT-qPCR change", ylab="MeRIP-seq L2FC",
          legend.title = '',legend='top',font.legend = c(8, "plain", "black")) 
dev.off()

pdf("fig4_merip_seq_vs_rtqpcr_all.pdf",height = 3.5,width=3.5)
seq.qpcr.r <- round(cor(merip.rt.qpcr.v.seq$Mean, merip.rt.qpcr.v.seq$MeRIPseq.change),2)
seq.qpcr.p <- formatC(cor.test(merip.rt.qpcr.v.seq$Mean, merip.rt.qpcr.v.seq$MeRIPseq.change)$p.value,format="e",digits=1)
ggplot(merip.rt.qpcr.v.seq,aes(x = Mean, y = MeRIPseq.change)) + geom_point(color="#db7889") + 
  geom_smooth(method = "lm",se=TRUE,colour="#8e1b64",fill="#8e5b7b",size=0.5) + 
  xlab("MeRIP-RT-qPCR L2FC") +
  ylab("MeRIP-seq L2FC") + scale_x_continuous(expand = c(0.01,0.01)) +
  theme_bw() + annotate(geom="text",x=-0.2, y=3,size=4,label=paste0('list(italic(R) == ',"`",seq.qpcr.r,"`",', italic(p) ==',"`",seq.qpcr.p,"`",")"), parse=TRUE) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()