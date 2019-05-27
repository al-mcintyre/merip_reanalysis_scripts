library(edgeR)
library(QNB)
library(DESeq2)
source('../fig2/metdiff_function_copy.R')

run.tool <- function(tool,cnts,meta,cond,full.results){
  if (tool == "DESeq2_GLM"){ 
    res <- run.deseq2(cnts,meta,l2fc = full.results)
  }else if (tool == "edgeR_GLM"){
      res <- run.edger(cnts,meta,l2fc = full.results)
  }else if (tool == "QNB"){
    run.qnb(cnts,meta,cond,l2fc = full.results)
  }else if (tool == "MeTDiff"){
    run.metdiff(cnts,meta,cond,l2fc = full.results)
    }
}

run.deseq2 <- function(cnts,meta,adjust=FALSE,l2fc=FALSE){
  inf.dds <- DESeqDataSetFromMatrix(countData = cnts,colData = meta,design = ~Condition+IP+Condition:IP)
  inf.dds.LRT <- DESeq(inf.dds,betaPrior=FALSE, test="LRT",
                       full=~Condition+IP+Condition:IP,reduced=~Condition+IP)    
  inf.dds.res <- results(inf.dds.LRT)
  results <- inf.dds.res$pvalue
  if (adjust){
    results <- inf.dds.res$padj
  }
  if (l2fc){
    results <- as.data.frame(cbind(inf.dds.res$log2FoldChange,inf.dds.res$pvalue,inf.dds.res$padj))
    colnames(results) <- c("l2fc","p","padj")
  }
  return(results)
}

run.qnb <- function(cnts,meta,label="Horner",adjust=FALSE,l2fc=FALSE){
  meth1 <- cnts[,which(meta$Condition == label & meta$IP == "IP")]
  meth2 <- cnts[,which(meta$Condition != label & meta$IP == "IP")]
  unmeth1 <- cnts[,which(meta$Condition == label & meta$IP == "input")]
  unmeth2 <- cnts[,which(meta$Condition != label & meta$IP == "input")]
  qnb.result.sim = qnbtest(meth1, meth2, unmeth1, unmeth2, mode="per-condition")
  results <- qnb.result.sim$pvalue
  if (adjust){
    results <- qnb.result.sim$padj
  }
  if (l2fc){
    results <- as.data.frame(cbind(qnb.result.sim$log2.RR,qnb.result.sim$log2.OR,qnb.result.sim$pvalue,qnb.result.sim$padj))
    colnames(results) <- c("l2rr","l2or","p","padj")
  }
  return(results)
}

run.metdiff <- function(cnts,meta,label="mockH",adjust=FALSE,l2fc=FALSE){
  meth1 <- cnts[,which(meta$Condition == label & meta$IP == "IP")]
  meth2 <- cnts[,which(meta$Condition != label & meta$IP == "IP")]
  unmeth1 <- cnts[,which(meta$Condition == label & meta$IP == "input")]
  unmeth2 <- cnts[,which(meta$Condition != label & meta$IP == "input")]
  metdiff.result <- diff.call.module(meth1,unmeth1,meth2,unmeth2)
  results <- metdiff.result$DIFF$pvalues
  if (adjust){
    results <- metdiff.result$DIFF$fdr
  }
  if (l2fc){
    results <- as.data.frame(cbind(metdiff.result$DIFF$fc,metdiff.result$DIFF$pvalues,metdiff.result$DIFF$fdr))
    colnames(results) <- c("fc","p","padj")
  }
  return(results)
}

run.edger <- function(cnts,meta,adjust=FALSE,l2fc=FALSE){
  #add count filter?
  er.design <- model.matrix(~meta$Condition+meta$IP+meta$Condition*meta$IP)
  er.dgelist <- DGEList(counts=cnts,group=meta$Condition) 
  er.dgelist <- estimateDisp(er.dgelist, design=er.design)
  er.fit <- glmFit(er.dgelist, er.design)
  er.lrt <- glmLRT(er.fit, coef=4)
  #hist(er.lrt$table$PValue)
  results <- er.lrt$table$PValue
  if (adjust){
    results <- p.adjust(er.lrt$table$PValue)
  }
  if (l2fc){
    results <- as.data.frame(cbind(er.lrt$table$logFC,er.lrt$table$PValue,p.adjust(er.lrt$table$PValue)))
    colnames(results) <- c("l2fc","p","padj")
  }
  return(results)
}

run.deseq2.multisample <- function(cnts,meta){
  inf.dds <- DESeqDataSetFromMatrix(countData = cnts,colData = meta,design = ~Tissue+Treatment+IP+Treatment:IP)
  inf.dds.LRT <- DESeq(inf.dds,betaPrior=FALSE, test="LRT",
                       full=~Tissue+Treatment+IP+Treatment:IP,reduced=~Tissue+Treatment+IP)    
  inf.dds.res <- results(inf.dds.LRT)
  return(inf.dds.res$pvalue) #adj
}

run.edger.multisample <- function(cnts,meta){
  er.design <- model.matrix(~meta$Tissue+meta$IP+meta$Treatment+meta$Treatment*meta$IP)
  er.dgelist <- DGEList(counts=cnts,group=meta$Treatment) 
  er.dgelist <- estimateDisp(er.dgelist, design=er.design)
  er.fit <- glmFit(er.dgelist, er.design)
  er.lrt <- glmLRT(er.fit, coef=5)
  return(er.lrt$table$PValue)
}

run.deseq.gene <- function(condition,control.cond,cnts,coldata,label){
  dds <- DESeqDataSetFromMatrix(cnts,coldata,formula(~Condition))
  dds$Condition <- factor(dds$Condition, levels=c(control.cond,condition))
  gene.col2check <- coldata$Condition
  dds$Condition <- droplevels(dds$Condition)
  gene.deseq <- DESeq(dds)
  gene.deseq <- results(gene.deseq)
  gene.results <- gene.deseq[,c("log2FoldChange","pvalue","padj")]
  colnames(gene.results) <- paste0(label,c("_l2fc","_p","_padj"))
  return(gene.results)
}