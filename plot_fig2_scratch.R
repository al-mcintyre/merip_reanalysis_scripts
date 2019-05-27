#expand GLMs to non-peak counts
exp.metadata.pk <- metadata
exp.metadata.pk$region <- "peak"
exp.metadata.npk <- metadata
exp.metadata.npk$region <- "nonpeak"
exp.metadata <- rbind(exp.metadata.pk,exp.metadata.npk)
gene.counts <- counts(macs2.peak.counts$gene.counts.nopeaks)
peak.genes <- macs2.peak.counts$peaks$exon

#for edgeR and DESeq2 log likelihood fits:   
deseq2.md <- get.deseq2.ests(subcounts,submeta)
deseq2.md <- deseq2.md[keep.peaks,]
plot(log10(deseq2.md[,"de.means"]),log10(deseq2.md[,"de.disps"]))

edger.md <- get.edgeR.ests(subcounts,submeta)
edger.md <- edger.md[keep.peaks,]
methods <- c("mle") #,"deseq2","edgeR")
mds <- data.frame("mle"=mle.md) #,"de"=deseq2.md,"er"=edger.md)

#3D plots once have subcounts with means and dispersions calculated
mle.var <- apply(subcounts,1,var)
mle.skew <- apply(subcounts,1,skewness)
plot(log10(mle.md[,"mu"]),log10(mle.var),xlab="log10(mean)",ylab="log10(var)")
plot(log10(mle.md[,"mu"]),mle.skew,xlab="log10(mean)",ylab="skewness")
plot(log10(mle.var),mle.skew,xlab="log10(var)",ylab="skewness")
scatterplot3d(x=log10(mle.md[,"mu"]), y=log10(mle.var), z=mle.skew, angle = 55)


#add error bars for downsampled data with permutations
neg.controls <- as.data.frame(cbind(c(mocks.2reps,mocks.3reps,mocks.full,DAA.2,DAA,shWTAP,METTL3kd,FTOoe),rep(c("DESeq2 GLM","edgeR GLM","MeTDiff","QNB"),8)))
neg.controls$sample <- c(rep("mocks (2 reps)",4),rep("mocks (3 reps)",4),rep("mocks (6 reps)",4),rep("DAA (2 reps)",4),rep("DAA (3 reps)",4),rep("shWTAP (2 reps)",4),rep("siMETTL3 (2 reps)",4),rep("FTOoe (2 reps)",4))
colnames(neg.controls) <- c("perc.p","Method","Sample")
neg.controls$perc.p <- as.numeric(as.character(neg.controls$perc.p))
summary <- summaryBy(perc.p ~ Method + Sample, neg.controls, FUN = c(mean, sd))
summary$perc.p.sd[which(is.na(summary$perc.p.sd))] <- 0
summary$Sample <- factor(summary$Sample,levels=c("mocks (2 reps)","mocks (3 reps)","mocks (6 reps)","DAA (2 reps)","DAA (3 reps)","shWTAP (2 reps)","siMETTL3 (2 reps)","FTOoe (2 reps)"))
pdf("controls_pvalues.pdf",width = 7, height=4)
ggplot(data=summary, aes(x=Method, y=perc.p.mean, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge()) + ylab("% peaks p < 0.05") +
  geom_errorbar(data = summary, mapping = aes(ymin = perc.p.mean - perc.p.sd, ymax = perc.p.mean + perc.p.sd),position=position_dodge(.9),width=0) +
  xlab("") +
  scale_fill_manual(values=c("#ff3c38","#ba324f","#6a0136","#d8f1a0","#00a878","#343e3d","#7a918d","#33673b")) + theme(panel.grid.major.y = element_line(colour = "lightgrey"), panel.background = element_blank())
dev.off()

perc.change.controls <- as.data.frame(cbind(c(mocks.perc.2reps,mocks.perc.3reps,mocks.perc.full,DAA.perc.2reps,DAA.perc.3reps,shWTAP.perc.2reps,siMETTL3.perc.2reps,FTOoe.perc.2reps),rep(c("DESeq2 GLM","edgeR GLM","MeTDiff","QNB"),8)))
perc.change.controls$sample <- c(rep("mocks (2 reps)",4),rep("mocks (3 reps)",4),rep("mocks (6 reps)",4),rep("DAA (2 reps)",4),rep("DAA (3 reps)",4),rep("shWTAP (2 reps)",4),rep("siMETTL3 (2 reps)",4),rep("FTOoe (2 reps)",4))
colnames(perc.change.controls) <- c("perc.change","Method","Sample")
perc.change.controls$perc.change <- as.numeric(as.character(perc.change.controls$perc.change))
summary <- summaryBy(perc.change ~ Method + Sample, perc.change.controls, FUN = c(mean, sd))
summary$perc.change.sd[which(is.na(summary$perc.change.sd))] <- 0
summary$Sample <- factor(summary$Sample,levels=c("mocks (2 reps)","mocks (3 reps)","mocks (6 reps)","DAA (2 reps)","DAA (3 reps)","shWTAP (2 reps)","siMETTL3 (2 reps)","FTOoe (2 reps)"))
pdf("controls_percent_changes.pdf",width = 7, height=4)
ggplot(data=summary, aes(x=Method, y=perc.change.mean, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge()) + ylab("% significant peaks\nwith change in %m6A <= -20") +
  geom_errorbar(data = summary, mapping = aes(ymin = perc.change.mean - perc.change.sd, ymax = perc.change.mean + perc.change.sd),position=position_dodge(.9),width=0) +
  xlab("") +
  scale_fill_manual(values=c("#ff3c38","#ba324f","#6a0136","#d8f1a0","#00a878","#343e3d","#7a918d","#33673b")) + theme(panel.grid.major.y = element_line(colour = "lightgrey"), panel.background = element_blank())
dev.off()

##Make simulated negative control sets by changing gene expression for mock data between labs
input.counts <- counts[,which(metadata$IP == "input")]
ip.counts <- counts[,which(metadata$IP == "IP")]
n <- dim(counts)[2]/2

#fit IP samples
mod = "Negative Binomial"
data <- "IP"
submeta <- metadata[which(metadata$IP == data),]
mle.md.ip <- apply(unname(ip.counts),1,fit.mle,model=mod)
mle.md.ip <- do.call(rbind, mle.md.ip)
keep.peaks <- which(complete.cases(mle.md.ip))
mle.md.ip <- mle.md.ip[keep.peaks,]
input.counts <- input.counts[keep.peaks,]
ip.counts <- ip.counts[keep.peaks,]

meandisp.fit.ip <- parametricDispersionFit(mle.md.ip[,1],1/mle.md.ip[,2])
plot(mle.md.ip[,1],1/mle.md.ip[,2],log = "xy",ylab="dispersion (MLE)",xlab="mean",pch=16)
pred <- meandisp.fit.ip(mle.md.ip[,1])
points(mle.md.ip[,1], pred, col='red',pch=16)

#fit input samples
data <- "input"
submeta <- metadata[which(metadata$IP == data),]
mle.md.input <- apply(unname(input.counts),1,fit.mle,model=mod)
mle.md.input <- do.call(rbind, mle.md.input)
keep.peaks <- which(complete.cases(mle.md.input))
mle.md.input <- mle.md.input[keep.peaks,]
input.counts <- input.counts[keep.peaks,]
ip.counts <- ip.counts[keep.peaks,]

meandisp.fit.input <- parametricDispersionFit(mle.md.input[,1],1/mle.md.input[,2])
plot(mle.md.input[,1],1/mle.md.input[,2],log = "xy",ylab="dispersion (MLE)",xlab="mean",pch=16)
pred <- meandisp.fit.input(mle.md.input[,1])
points(mle.md.input[,1], pred, col='red',pch=16)

iterations <- 3
intervention <- "unsorted" #vs. sorted, 2x, 5x, 10x
factor <- 2

for (j in 2:iterations){
  mod.counts <- cbind(ip.counts,input.counts) #matrix(ncol=dim(counts)[2], nrow=dim(counts)[1])
  mod.inds <- sample(1:dim(ip.counts)[1],1000)
  for (i in mod.inds){
    print(i)
    row <- input.counts[i,]
    ip.row <- ip.counts[i,]
    if (intervention == "sorted"){ #separate highest and lowest mock samples by reordering based on inputs
      ord <- order(row)
      if (i%%2 == 0){
        ord <- rev(ord)
      }
      mod.counts[i,] <- c(ip.row[ord],row[ord])
    } else if (intervention == "2x"){ #multiply mean by factor and calculate new dispersion based on model fit with all means/dispersions (separately for IPs and inputs)
      md.1.ip <- mle.md.ip[i,]
      md.2.ip <- c(md.1[1]*factor, 1/meandisp.fit.ip(md.1.ip[1]*factor))
      md.1.input <- mle.md.input[i,]
      md.2.input <- c(md.1.input[1]*factor, 1/meandisp.fit.input(md.1.input[1]*factor))
      #could also pick more typical dispersion for the data with unmodified mean
      mod.counts[i,] <- c(sim.data(md.1.ip,6,"Negative Binomial"),sim.data(md.2.ip,6,"Negative Binomial"),sim.data(md.1.input,6,"Negative Binomial"),sim.data(md.2.input,6,"Negative Binomial"))
    } 
    #else if (intervention == "unsorted"){
    #}
  }
  deseq.p.mod <- run.deseq2(mod.counts,metadata)
  edger.p.mod <- run.edger(mod.counts,metadata)
  qnb.p.mod <- run.qnb(mod.counts,metadata,"mockH")
  metdiff.p.mod <- run.metdiff(mod.counts,metadata,"mockH")
  mocks.mod <- c(length(which(deseq.p.mod[mod.inds] < 0.05)),length(which(edger.p.mod[mod.inds] < 0.05)),length(which(metdiff.p.mod[mod.inds] < 0.05)),length(which(qnb.p.mod[mod.inds] < 0.05)))*100/length(deseq.p.mod[mod.inds])
  description.mod <- c("","","")
  pvals <- as.data.frame(cbind(mocks.mod,c("DESeq2 GLM","edgeR GLM","MeTDiff","QNB"),rep(paste(intervention,"(6 reps)"),4)))
  colnames(pvals) <- c("perc.p","Method","Sample")
  if (exists("mod.pval")){
    mod.pval <- rbind(mod.pval,pvals)
  }else{
    mod.pval <- pvals
  }  
}

#add error bars for downsampled data with permutations
#colnames(mod.pval) <- c("perc.p","Method","Sample")
mod.pval$perc.p <- as.numeric(as.character(mod.pval$perc.p))
#mod.pval <- mod.pval[9:dim(mod.pval)[1],]
summary <- summaryBy(perc.p ~ Method + Sample, mod.pval, FUN = c(mean, sd))
summary$perc.p.sd[which(is.na(summary$perc.p.sd))] <- 0
summary$Sample <- factor(summary$Sample,levels=c("unsorted (6 reps)","sorted (6 reps)"))
pdf("mocks_pvalues.pdf",width = 5.3, height=4)
ggplot(data=summary, aes(x=Method, y=perc.p.mean, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge()) + ylab("% peaks p < 0.05") +
  geom_errorbar(data = summary, mapping = aes(ymin = perc.p.mean - perc.p.sd, ymax = perc.p.mean + perc.p.sd),position=position_dodge(.9),width=0) +
  ylim(c(0,25)) + xlab("") +
  scale_fill_manual(values=c("#6a0136","#5b8e7d","#f4a259","#8cb369","#f4e285")) + theme(panel.grid.major.y = element_line(colour = "lightgrey"), panel.background = element_blank())
dev.off()

#redo with only three replicates per condition

#could also try making mock positive control by pairing highest IP counts with lowest input counts and vice versa

#run DESeq2, edgeR, MeTDiff, QNB on the mock data

#borrowing disperision fit function from DESeq2
parametricDispersionFit <- function( means, disps ) {
  coefs <- c( .1, 1 )
  iter <- 0
  while(TRUE) {
    residuals <- disps / ( coefs[1] + coefs[2] / means )
    good <- which( (residuals > 1e-4) & (residuals < 15) )
    # check for glm convergence below to exit while-loop
    suppressWarnings({fit <- glm( disps[good] ~ I(1/means[good]),
                                  family=Gamma(link="identity"), start=coefs )})
    oldcoefs <- coefs
    coefs <- coefficients(fit)
    if ( !all( coefs > 0 ) )
      stop(simpleError("parametric dispersion fit failed"))
    if ( ( sum( log( coefs / oldcoefs )^2 ) < 1e-6 )  & fit$converged )
      break
    iter <- iter + 1
    if ( iter > 10 ) 
      stop(simpleError("dispersion fit did not converge"))
  }
  names( coefs ) <- c( "asymptDisp", "extraPois" )
  ans <- function(q) coefs[1] + coefs[2] / q
  attr( ans, "coefficients" ) <- coefs
  ans
}
