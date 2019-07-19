get.gtf <- function(species){
  if (species == "hg38"){ return("") #fill in hg38 gtf!
  } else if (species == "mm10"){
    return("") #fill in mm10 gtf!
  }
}

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


count.matches <- function(pat, vec) sapply(regmatches(vec, gregexpr(pat, vec,perl=TRUE)), length)

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
