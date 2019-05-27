.help.digamma <- function(xx,alpha){
  Tm <- dim(xx)
  TT <- Tm[1]
  m <- Tm[2]
  
  res <- matrix(0,TT,1)
  
  for (ii in 1:m){
    res <- res + digamma(xx[,ii] + alpha)
  }
  return(res)
}

.help.trigamma <- function(xx,alpha){
  Tm <- dim(xx)
  if (is.null(Tm)){
    TT <- length(xx)
    m <- 1
  }else{
    TT <- Tm[1]
    m <- Tm[2]
  }
  res <- matrix(0,TT,1)
  
  for (ii in 1:m){
    res <- res + trigamma(xx[,ii] + alpha)
  }
  return(res)
  
}

.help.postprob <- function(dxx,dyy,dnn,xx,yy,alpha,beta){
  
  N <- length(alpha)
  Tm <- dim(xx)
  if (is.null(Tm)){
    TT <- length(xx)
    m <- 1
  }else{
    TT <- Tm[1]
    m <- Tm[2]
  }
  res <- matrix(0,TT,N)
  
  for (ii in 1:m){
    dnx <- as.matrix(dxx[,ii])
    dny <- as.matrix(dyy[,ii])
    dn <- as.matrix(dnn[,ii])
    x <- as.matrix(xx[,ii])
    y <- as.matrix(yy[,ii])
    res <- res + (dn-dnx-dny) %*% matrix(1,1,N) + lgamma(x %*% matrix(1,1,N) + matrix(1,TT) %*% alpha) - 
      lgamma(matrix(1,TT) %*% (alpha+beta) + (x+y) %*% matrix(1,1,N)) +
      lgamma(y %*% matrix(1,1,N) + matrix(1,TT) %*% beta) + lgamma(matrix(1,TT) %*% (alpha+beta)) - 
      lgamma(matrix(1,TT) %*% alpha) - lgamma(matrix(1,TT) %*% beta)
  }
  res <- exp(res)
}

.help.factorial <- function(count){
  #compute the log(count!)
  cm = max(count)
  if (is.null(ncol(count))){
    D <- 1
  }else{
    D <- ncol(count)
  }
  
  if(cm > 50000){
    dnorm <- as.matrix(lgamma(data.matrix(count+1)))
  }
  else{
    tmp  <- cumsum(rbind(0,log(as.matrix(1:max(count)))))
    dnorm <- matrix(tmp[data.matrix(count+1)],ncol=D)
  }
}

.betabinomial.lh <- function(x,y,Nit=40,Npara=1e-9){
  #   x <- as.matrix(x[peak,])
  #   y <- as.matrix(y[peak,])
  N <- 2 # number of states
  J <- matrix(0,N,1)
  H <- matrix(0,N,N)
  T <- nrow(x)
  
  IP_mean <- rowMeans(x)
  INPUT_mean <- rowMeans(y)
  nip = ncol(x)
  nin = ncol(y)
  # if the dimension for x and y does not match
  if (nip > nin) {
    avg_input <- round(matrix(rep(INPUT_mean,nip-nin),ncol=nip-nin))
    y <- cbind(y,avg_input) 
  }
  else if (nip < nin){
    avg_ip <- matrix(rep(IP_mean,nin-nip),ncol=nin-nip)
    x <- cbind(x,avg_ip) 
  }
  
  
  
  n <- x + y
  m <- ncol(x)
  rr <- x/n
  # parameters initialization
  #   exx <- colMeans(rr)
  #   vxx <- apply(rr,2,sd)^2
  #   alpha <- mean( exx*(exx*(1-exx)/vxx-1) )
  #   beta <- mean( (1-exx)*(exx*(1-exx)/vxx-1) )
  
  # use another method to initialize
  p1_e <- exp(sum( log(rr) )/(T*m))
  p2_e <- exp(sum( log(1-rr)/(T*m) ))
  alpha <- 1/2 *(1-p2_e)/(1-p1_e-p2_e ) # to avoid 0
  beta <-  1/2 *(1-p1_e)/(1-p1_e-p2_e )
  c = rbind(alpha,beta)
  
  # add break condition to avoid alpha is na   
  if ( !any(is.finite(beta)) | !is.finite(alpha) | any(beta <= 0) | any(alpha<= 0) ){
    return(list(logl=rnorm(1)*10000,alpha=c(1,1),beta=c(1,1))) 
  }
  
  
  
  for (nit in 1:Nit){
    J[1] <- T*digamma(sum(c))*m - sum( .help.digamma(as.matrix(n),sum(c)) ) + sum( .help.digamma(as.matrix(x),c[1]) ) - T*digamma(c[1])*m
    J[2] <- T*digamma(sum(c))*m - sum( .help.digamma(as.matrix(n),sum(c)) ) + sum( .help.digamma(as.matrix(y),c[2]) ) - T*digamma(c[2])*m
    
    H[1,1] <- T*trigamma(sum(c))*m - sum(.help.trigamma(as.matrix(n),sum(c))) + sum(.help.trigamma(as.matrix(x),c[1])) - T*trigamma(c[1])*m
    H[2,2] <- T*trigamma(sum(c))*m - sum(.help.trigamma(as.matrix(n),sum(c))) + sum(.help.trigamma(as.matrix(y),c[2]))  - T*trigamma(c[2])*m    
    H[1,2] <- T*trigamma(sum(c))*m - sum(.help.trigamma(as.matrix(n),sum(c)))
    H[2,1] <- H[1,2]
    eigvalue <- eigen(H)$values
    
    if ( (any(beta < Npara)) | (any(alpha < Npara)) 
         | abs(eigvalue[1]/eigvalue[2]) > 1e12 | abs(eigvalue[1]/eigvalue[2]) < 1e-12
         | any(eigvalue==0) ){   break  }
    
    #     tmp_step <- -solve(H,tol=1e-20) %*% J
    tmp_step <- -solve(H, J) # using newton smoothing
    tmp <- c + tmp_step
    while(any(tmp <= 0)){
      #       warning(sprintf("Could not update the Newton step ...\n"))
      tmp_step <- tmp_step / 20
      tmp <- c + tmp_step
    }
    c <- tmp
    
  }
  #   caculate the likelihood
  alpha <- c[1]
  beta <- c[2]
  dnx <- .help.factorial(x)
  dny <- .help.factorial(y)
  dn <- .help.factorial(n)
  prob <- .help.postprob(dnx,dny,dn,x,y,alpha,beta)
  return(list(logl=sum(log(prob)),alpha=alpha,beta=beta))
  
}

# merge and compare two conditions
diff.call.module <- function(meth1,unmeth1,meth2,unmeth2){
  #x = untreated IP, y = untreated input, xx = treated IP, yy = treated input
  no_peak=length(meth1[,1]) #PEAK$loci2peak_merged[,1])
  pvalues <- rep(1,no_peak)
  log.fc <- rep(0,no_peak)
  for (ipeak in 1:no_peak) {
    if (ipeak%%1000 == 0){print(ipeak)}
    x = t(as.array(meth1[ipeak,]))
    y = t(as.matrix(unmeth1[ipeak,]))
    xx = t(as.matrix(meth2[ipeak,]))
    yy = t(as.matrix(unmeth2[ipeak,]))
    xxx = cbind(x,xx)
    yyy = cbind(y,yy)
    #BBtest
    logl1 <- .betabinomial.lh(x,y+1)
    logl2 <- .betabinomial.lh(xx,yy+1)
    logl3 <- .betabinomial.lh(xxx,yyy+1)
    tst <- (logl1$logl+logl2$logl-logl3$logl)*2
    pvalues[ipeak] <- 1 - pchisq(tst,2)
    log.fc[ipeak] <- log2( (sum(xx)+1)/(1+sum(yy)) * (1+sum(y))/(1+sum(x)) ) 
    
  }
  
  p <- pvalues
  fdr <- p.adjust(pvalues,method='fdr')
  #log.fdr <- pmax(log.fdr,-1000)
  #log.p <- pmax(log.p,-1000)
  DIFF <- list(fdr=fdr,pvalues=p,fc=log.fc)
  # result
  result =list()
  result$DIFF = DIFF
  #   result$peak_reads_count = list(untreated_ip=untreated_ip,
  #                                  untreated_input=untreated_input,
  #                                  treated_ip=treated_ip,
  #                                  treated_input=treated_input,
  #                                  untreated_ip_total=untreated_ip_total,
  #                                  untreated_input_total=untreated_input_total,
  #                                  treated_ip_total=treated_ip_total,
  #                                  treated_input_total=treated_input_total,
  #                                  sample_peak_reads_count=peak_reads_count,
  #                                  sample_total_reads_count = sample_total,
  #                                  sample_id=SAMPLE_ID)
  return(result)
  
}
