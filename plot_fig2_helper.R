library(CLIPanalyze)

calc.ll <- function(c.md,model="Negative Binomial"){
  if (model=="Negative Binomial"){
    len <- length(c.md)
    m <- c.md[len-1] #rep(c.md[len-1],len-2)
    d <- c.md[len] #rep(c.md[len],len-2)
    k <- c.md[1:(len-2)]
    log.likelihood <- sum(dnbinom(x=k,size=d,mu=m,log=TRUE)) #sum(mapply(calc.ll,k,m,d)) OR IS SIZE=1/d??
    return(log.likelihood)
  } else if (model=="Poisson"){
    len <- length(c.md)
    l <- c.md[len]
    k <- c.md[1:(len-1)]
    log.likelihood <- sum(dpois(x=k,lambda=l,log=TRUE)) 
    return(log.likelihood)    
  }
}

annot.strand <- function(peaks.etc,anno=annot){
  keys <- keys(anno$txdb, keytype="GENEID")
  cols <- c("EXONSTRAND")
  txdb <- anno$txdb
  geneid.strand <- select(txdb, keys=keys, columns=cols, keytype="GENEID") #select from GenomicFeatures, cannot have dplyr loaded (in MeTDiff?)
  gene.df <- as.data.frame(anno$genenames)
  row.names(gene.df) <- gene.df$gene_id
  row.names(geneid.strand) <- make.names(gene.df[geneid.strand$GENEID,"gene_name"],unique = TRUE)
  
  overlaps <- as.data.frame(peaks.etc$overlaps)
  row.names(overlaps) <- make.names(paste0("peak",overlaps$name),unique = TRUE)
  peak.main.genes <- overlaps[names(peaks.etc$peaks),"gene_name"]
  
  new.strands <- geneid.strand[peak.main.genes,"EXONSTRAND"]
  new.strands[is.na(new.strands)] <- "*"

  alt.geneid <- geneid.strand[peak.main.genes,"GENEID"]
  added.annotations <- cbind(new.strands,peak.main.genes,alt.geneid)
  return(added.annotations)
}

clip.annotate <- function(m2peaks,anno,annot.order){
  names(m2peaks) <- c(paste0("peak",c(1:length(m2peaks))))
  m2peaks.anno <- annotateFeatures(m2peaks, txdb = anno$txdb,
                                   genenames = anno$genenames,
                                   mirnas = anno$mirnas,
                                   utr5.extend = 2000,
                                   utr3.extend = 2000)
  names(m2peaks.anno) <- c("peaks", "overlaps")
  m2peaks.anno$peaks <- prioritizeAnnotations(m2peaks.anno$peaks,feature.order = annot.order) 
  g.s <- annot.strand(m2peaks.anno)
  strand(m2peaks.anno$peaks) <- g.s[,1]
  m2peaks.anno$peaks$maingene <- g.s[,2]
  m2peaks.anno$peaks$geneid <- g.s[,3]
  return(m2peaks.anno)
}

import.peaks <- function(files,anno=annot,annot.order=annotation.order){
  #remove(m2peaks)
  for (file in files){
    if (exists("m2peaks")){
      m2peaks <- union(m2peaks,import(file,format="BED"))
    }else{
      m2peaks <- import(file,format="BED")
    }
  }
  m2peaks <- reduce(m2peaks)
  m2peaks <- keepStandardChromosomes(m2peaks,pruning.mode = "coarse")
  #annotate macs2 peaks
  m2peaks.anno <- clip.annotate(m2peaks,anno,annot.order)
  remove(m2peaks)
  return(m2peaks.anno)
}

get.counts <- function(meta,m2peaks,bsg=bsgenome){
  bamlist <- c(as.vector(meta$bam))
  sample.names <-meta$SampleID
  m2.peak.counts <- countReads(m2peaks$peaks,bamlist,sample.names,genomeTag,stranded=FALSE,exons.only=TRUE,paired.end = TRUE) #consider paired end adjustment?
  strand(m2.peak.counts$peaks)[strand(m2.peak.counts$peaks) == "*"] <- "+" #necessary - can't mix * with +/- when getting sequences
  m2.peak.counts$peaks$seqs <- get.seqs(bsg,m2.peak.counts$peaks)
  #count % peaks containing DRACN motifs
  m2.peak.counts$peaks$dracn <- grepl('[AGT][AG]AC[ACGT]',m2.peak.counts$peaks$seqs)
  m2.peak.counts$peaks$dracn[m2.peak.counts$peaks$dracn == TRUE] <- "DRACN"
  m2.peak.counts$peaks$dracn[m2.peak.counts$peaks$dracn == FALSE] <- "non-motif"
  m2.peak.counts$peaks$length <- width(m2.peak.counts$peaks) 
  paste("% DRACN peaks",length(grep('[AGT][AG]AC[ACGT]',m2.peak.counts$peaks$seqs))*100/length(m2.peak.counts$peaks$seqs))
  return(m2.peak.counts)
}

count.matches <- function(pat, vec) sapply(regmatches(vec, gregexpr(pat, vec,perl=TRUE)), length)

finish.annotation <- function(meta,m2.peaks,bsg=bsgenome){
  sample.names <-meta$SampleID
  strand(m2.peaks$peaks)[strand(m2.peaks$peaks) == "*"] <- "+" #necessary - can't mix * with +/- when getting sequences
  m2.peaks$peaks$seqs <- get.seqs(bsg,m2.peaks$peaks)
  #count % peaks containing DRACN motifs
  m2.peaks$peaks$dracn <- count.matches("[AGT][AG]AC(?=[ACGT])",m2.peaks$peaks$seqs) #grepl('[AGT][AG]AC[ACGT]',m2.peaks$peaks$seqs)
  #m2.peaks$peaks$dracn[m2.peaks$peaks$dracn == TRUE] <- "DRACN"
  #m2.peaks$peaks$dracn[m2.peaks$peaks$dracn == FALSE] <- "non-motif"
  m2.peaks$peaks$ggacu <- count.matches("GGACT",m2.peaks$peaks$seqs) #grepl('[AGT][AG]AC[ACGT]',m2.peaks$peaks$seqs)
  #m2.peaks$peaks$ggacu <- grepl('GGACT',m2.peaks$peaks$seqs)
  #m2.peaks$peaks$ggacu[m2.peaks$peaks$ggacu == TRUE] <- "GGACU"
  #m2.peaks$peaks$ggacu[m2.peaks$peaks$ggacu == FALSE] <- "non-ggacu"
  m2.peaks$peaks$length <- width(m2.peaks$peaks) 
  paste("% DRACN peaks",length(grep('[AGT][AG]AC[ACGT]',m2.peaks$peaks$seqs))*100/length(m2.peaks$peaks$seqs))
  return(m2.peaks)
}

get.gene.counts <- function(meta,gtffi){
  bamlist <- c(as.vector(meta$bam)) 
  sample.names <-meta$SampleID
  gene.counts <- Rsubread::featureCounts(files = bamlist, GTF.featureType="exon",GTF.attrType="gene_id",
                         annot.ext=gtffi, isGTFAnnotationFile=TRUE,isPairedEnd = TRUE )
  colnames(gene.counts$counts) <- sample.names
  return(gene.counts$counts)
}

#build DESeq2 NB model
get.deseq2.ests <- function(subcounts,submeta){
  inf.dds <- DESeqDataSetFromMatrix(countData = subcounts,colData = submeta,design = ~1) 
  inf.dds <- estimateSizeFactors(inf.dds)
  inf.dds <- estimateDispersions(inf.dds)
  #hist(inf.dds.res$pvalue,col="#66D7D1",xlab="GLM p values",main="")
  de.disps <- 1/dispersions(inf.dds)
  #sfs <- sizeFactors(inf.dds.LRT)
  de.means <- rowData(inf.dds)$baseMean
  #plotDispEsts(inf.dds.LRT)
  return(cbind(de.means,de.disps))
}

#build edgeR NB model
get.edgeR.ests <- function(subcounts,submeta){
  #er.design <- model.matrix(~rep(1,dim(subcounts)[2])) #metadata$Lab+metadata$IP+metadata$Lab*metadata$IP)
  er.dgelist <- DGEList(counts=counts) #,group=metadata$Lab)
  er.dgelist <- estimateDisp(er.dgelist) #,er.design #this returns tagwise dispersions squeezed towards the common dispersion
  er.disps <- 1/er.dgelist$tagwise.dispersion
  er.means <- rowMeans(subcounts)
  #er.out <- addPriorCount(counts)
  #offsets <- er.out$offset[1,] #offsets correspond to log-transformed modified library sizes
  #mean(adjustedProfileLik(1/er.disps,counts,er.design,offset=NULL,adjust = FALSE)) #no difference with and without offsets, only works with non-intercept design
  return(cbind(er.means,er.disps))
}

#refit negative binomial or poisson model
fit.mle <- function(k,model="Poisson"){
  n <- length(k)
  tryCatch({
    fd <- MASS::fitdistr(k,model)
    if (model == "Negative Binomial"){
      return(c(fd$estimate['mu'],fd$estimate['size']))
    } else if (model == "Poisson"){
      return(fd$estimate['lambda'])
    }
  }, error = function(err){ 
    return(NA) #rep(0,n) should remove these peaks from both sample and simulated data
  })
}

#simulate data
sim.data <- function(md,n = nexps,model="Poisson"){
  if (model == "Negative Binomial"){
    return(rnbinom(n=n,size=as.numeric(md[2]),mu=as.numeric(md[1])))
  } else if (model == "Poisson"){
    return(rpois(n=n,lambda=md))
  }
}

perc.m6A.down <- function(cnts,pvalues){
  num.samples <- dim(cnts)[2]
  sig.set <- which(pvalues < 0.05)
  tot.below.0.05 <- length(sig.set)
  inp.cnts <- cnts[sig.set,((num.samples/2)+1):num.samples]
  ip.cnts <- cnts[sig.set,1:(num.samples/2)]
  perc.m6A <- 100*ip.cnts/(ip.cnts+inp.cnts)
  #assumes modified condition always precedes control
  delta.perc.m6A <- rowMeans(perc.m6A[,1:(num.samples/4)]) - rowMeans(perc.m6A[,((num.samples/4)+1):(num.samples/2)]) 
  return(length(which(delta.perc.m6A <= -20))*100/tot.below.0.05)
}