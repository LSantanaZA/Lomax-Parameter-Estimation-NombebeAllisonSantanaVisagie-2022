rm(list = ls())
cat("\014")
graphics.off()

set.seed(123)

#----------------------------
#Install and load packages:
#----------------------------
instll <- c("xtable", "maxLik", "Renext","miscTools","plm","EnvStats","MASS","actuar","evmix")
instlld <- dimnames(installed.packages())[[1]]
for (i in 1:length(instll)) {
  if (!(instll[i] %in% instlld)) {
    install.packages(instll[i])
  }
  library(instll[i], character.only = TRUE)
}

#----------------------------
#raw moments
muk <- function(x, k) {
  mean(x ^ k)
}

#----------------------------
#Lomax distribution function
pLomax <- function(q, beta, sigma) {
  1 - (sigma / (q + sigma)) ^ beta
}
#----------------------------
#Lomax quantile function
qLomax <- function(x, beta, sigma) {
  sigma * ((1 - x) ^ (-1 / beta) - 1)
}
#----------------------------
#Lomax density function
dLomax <- function(x, beta, sigma) {
  (beta * sigma ^ beta) / ((x + sigma) ^ (beta + 1))
}
#----------------------------
#Lomax random number generator
rLomax <- function(n, beta, sigma) {
  U <- runif(n)
  qLomax(U, beta, sigma)
}
#----------------------------
#Log-likelihood function (Lomax) # used in maxLik()
llik <- function(param, x) {
  beta <- param[1]
  sigma <- param[2]
  if (beta > 0 & sigma > 0) {
    ans <- log(beta) - log(sigma) - (1 + beta) * (log(1 + x / sigma))
  } else{
    ans <- NA
  }
}

#----------------------------
#Different phi function choices
phi.kl <- function(t) {
  -log(t)
} #Should be t*log(t), but there are numerical complications
#----------------------------
phi.chisq <- function(t) {
  (1 - t) ^ 2
}
#----------------------------
phi.tv <- function(t) {
  abs(t - 1)
}
#--------------------------------------------------
#Method of moments estimator
MME <- function(x) {
  n <- length(x)
  x <- sort(x)
  mu1 <- mean(x)
  mu2 <- mean(x^2)
  b <- 2*(mu2 - mu1^2)/(mu2 - 2* mu1 ^ 2)
  s <- mu1 * (b - 1)
  ans2 <- c(b, s)
  ans3 <- list(beta = ans2[1], sigma = ans2[2])
  return(ans3)
}
#----------------------------
#MME beta estimator (separate)
MMEbeta <- function(x) {
  n <- length(x)
  x <- sort(x)
  mu1 <- mean(x)
  mu2 <- mean(x^2)
  b <- 2*(mu2 - mu1^2)/(mu2 - 2*mu1^2)
  #b <- 2*(mu1 - mu1^2)/(mu2 - 2*mu2^2)
  #b <- 2 * (mu2 - mu1 ^ 2) / (mu2 - 2 * mu1 ^ 2)
  ans2 <- b
  return(ans2)
}
#----------------------------
#MME sigma estimate (separate)
MMEsigma <- function(x) {
  n <- length(x)
  x <- sort(x)
  mu1 <- mean(x)
  mu2 <- mean(x^2)
  b <- 2*(mu2 - mu1^2)/(mu2 - 2*mu1^2)
  s <- mu1*(b-1)
  ans2 <- s
  return(ans2)
}
#----------------------------
#MLE.b: Bias correction Method of moments estimator
MME.b <- function(x, B = 100) {
  n <- length(x)
  x <- sort(x)
  bstar <- numeric(B)
  sstar <- numeric(B)
  P <- matrix(0,ncol=n,nrow=B)
  bb <- 1
  while(bb<=B){
    j <- sample(1:n,n,replace=TRUE)
    Xstar <- x[j]
    CV.flag <-  (sqrt((n - 1)/n)*sd(Xstar)/mean(Xstar) < 1) # Check for existence of estimator
    if(!CV.flag){
      bstar[bb] <- MMEbeta(Xstar)
      sstar[bb] <- MMEsigma(Xstar)
      j <- factor(j, levels=1:n)
      P[bb,]<-table(j)/n
      bb<- bb + 1
    }
  }
  Pbar <- apply(P,2,mean)
  mutilde1 <- sum(Pbar*x)
  mutilde2 <- sum(Pbar*x^2)
  btilde <- 2*(mutilde2 - mutilde1^2)/(mutilde2 - 2*mutilde1^2)
  stilde <- mutilde1*(btilde-1)
  scor <-MMEsigma(x) - mean(sstar) + stilde
  bcor <-MMEbeta(x)  - mean(bstar) + btilde
  ans3 <- list(beta = bcor, sigma = scor)
  return(ans3)
}


#----------------------------
#L-moments estimator:
#--------------------------------------------------
#Hosking, J. R. M. (1990). L-moments: Analysis and estimation of distributions using linear combinations of order statistics. Journal of the Royal Statistical Society. Series B (Methodological), 52, 105-124.
#----------------------------
l2 = function(X) {
  n <- length(X)
  X <- sort(X)
  nchoose2 <- choose(n, 2)
  j <- sequence(1:n)
  i <- rep(1:n, 1:n)
  Answer <- sum(X[i] - X[j]) / (2 * nchoose2)
  return(Answer)
}
#LME: The individual \hat{\beta} function
betahat.L <- function(x) {
  L2 <- l2(x)
  L1 <- mean(x)
  answer <- L2 / (2 * L2 - L1)
  return(answer)
}
#----------------------------
#LME: The individual \hat{\sigma} function
sigmahat.L <- function(x) {
  L2 <- l2(x)
  L1 <- mean(x)
  answer <- (L1 ^ 2 - L1 * L2) / (2 * L2 - L1)
  return(answer)
}
#----------------------------
#The L-moment estimator function:
LME <- function(x) {
  L2 <- l2(x)
  L1 <- mean(x)
  bhat.L <- L2 / (2 * L2 - L1)
  shat.L <- (L1 ^ 2 - L1 * L2) / (2 * L2 - L1)
  ans3 <- list(beta = bhat.L, sigma = shat.L)
  return(ans3)
}
#--------------------------------------------------
#Minimum Distance Estimator (CvM)
MDE.CvM <- function(x) {
  n <- length(x)
  x <- sort(x)
  MD.CvM <- function(param) {
    beta <- param[1]
    sigma <- param[2]
    ans1 <-   1/(12*n) + sum((pLomax(x,beta,sigma)-(2*(1:n)-1)/(2*n))^2)
    return(ans1)
  }
  lme <- LME(x)
  b <- lme$beta
  s <- lme$sigma
  ans2 <- optim(c(b, s), MD.CvM, method = "BFGS")
  ans3 <- list(beta=ans2$par[1],sigma=ans2$par[2],convergence=ans2$convergence)
  return(ans3)
}
#----------------------------
#MDE.SD/MDE.LS: Minimum Distance Estimator (LS or SD in the paper)
MDE.LS <- function(x) {
  n <- length(x)
  x <- sort(x)
  MD.LS <- function(param) {
    beta <- param[1]
    sigma <- param[2]
    ans2 <- (1 / n) * sum((pLomax(x, beta, sigma) - (1:n) / (n + 1)) ^ 2)
    return(ans2)
  }
  lme <- LME(x)
  b <- lme$beta
  s <- lme$sigma
  ans2 <- optim(c(b, s), MD.LS, method = "BFGS")
  ans3 <- list(beta=ans2$par[1],sigma=ans2$par[2],convergence=ans2$convergence)
  return(ans3)
}
#----------------------------
#MDE.KL: Minimum Distance Estimator (Phi=KL)
MDE.Phi.kl <- function(x) {
  n <- length(x)
  x <- sort(x)
  MD.Phi.kl <- function(param) {
    beta <- param[1]
    sigma <- param[2]
    fh <- density(x)
    fhat <- approxfun(fh$x, fh$y)
    fdg <- dLomax(x, beta, sigma) / fhat(x)
    mean(phi.kl(fdg))
  }
  lme <- LME(x)
  b <- lme$beta
  s <- lme$sigma
  ans2 <- optim(c(b, s), MD.Phi.kl, method = "Nelder-Mead")
  ans3 <- list(beta = ans2$par[1], sigma = ans2$par[2], convergence=ans2$convergence)
  return(ans3)
}
#----------------------------
#MDE.chisq: Minimum Distance Estimator (Phi=ChiSq)
MDE.Phi.chisq <- function(x) {
  n <- length(x)
  x <- sort(x)
  MD.Phi.chisq <- function(param) {
    beta <- param[1]
    sigma <- param[2]
    fh <- density(x)
    fhat <- approxfun(fh$x, fh$y)
    fdg <- dLomax(x, beta, sigma) / fhat(x)
    gdf <- fdg ^ (-1)
    mean(gdf * phi.chisq(fdg))
  }
  lme <- LME(x)
  b <- lme$beta
  s <- lme$sigma
  ans2 <- optim(c(b, s), MD.Phi.chisq, method = "Nelder-Mead")
  ans3 <- list(beta = ans2$par[1], sigma = ans2$par[2], convergence=ans2$convergence)
  return(ans3)
}
#----------------------------
#MDE.TV: Minimum Distance Estimator (Phi=Total Variation)
MDE.Phi.tv <- function(x) {
  n <- length(x)
  x <- sort(x)
  MD.Phi.tv <- function(param) {
    beta <- param[1]
    sigma <- param[2]
    fh <- density(x)
    fhat <- approxfun(fh$x, fh$y)
    fdg <- dLomax(x, beta, sigma) / fhat(x)
    gdf <- fdg ^ (-1)
    mean(gdf * phi.tv(fdg))
  }
  lme <- LME(x)
  b <- lme$beta
  s <- lme$sigma
  ans2 <- optim(c(b, s), MD.Phi.tv, method = "Nelder-Mead")
  ans3 <- list(beta = ans2$par[1], sigma = ans2$par[2], convergence=ans2$convergence)
  return(ans3)
}
#----------------------------------------------------
#MLE: Maximum Likelihood Estimator
MLE <- function(x){
  ans <- flomax(x)
  ans3 <- list(beta = as.numeric(ans[[1]][1]), sigma = as.numeric(ans[[1]][2]), returnCode = as.numeric(!ans$cvg))
  return(ans3)
}
#----------------------------
#Bias Adjusted Maximum Likelihood Estimator
MLE.b <- function(x, ruleofthumb = TRUE) {
  n <- length(x)
  x <- sort(x)
  mle <- MLE(x)
  b <- mle$beta
  s <- mle$sigma
  #D <- as.numeric(151 <= n & n <= 500)
  D <- as.numeric(151 <= n)
  bstar <-
    -0.991404 + 0.027773 * n - 0.0000328 * n ^ 2 + 0.64541 * log(n) - 3.617656 *
    D - 0.014898 * D * n + 0.0000293 * D * n ^ 2 + 1.036946 * D * log(n)
  if (b < bstar & ruleofthumb == TRUE) {
    K <- n*matrix(c(b/((b + 2)*s^2),
                    -1/((b + 1)*s),
                    -1/((b+1) * s),
                    1/b^2), ncol = 2)
    A <- n*matrix(c(2*b/((b+2)*(b+3)*s^3),
                    -1/((b+1)*(b+2)*s^2),
                    b/((s^2)*((b+2)^2)),
                    -1/(s*(b+1)^2),
                    -1/((s^2)*(b+1)*(b+2)),
                    0,
                    -1/(s*(b+1)^2),
                    1/b^3),
                  byrow = TRUE,
                  ncol = 4)
    ans <-c(mle$sigma, mle$beta)-as.vector(solve(K)%*%A%*%as.vector(solve(K)))
  } else{
    ans <- c(mle$sigma, mle$beta)
  }
  return(list(beta = ans[2], sigma = ans[1]))
}
#----------------------------
#PWM: functions for the PWME
M1u0 = function(X, u) {
  #u can ONLY be 0 or 1
  n <- length(X)
  X <- sort(X)
  i <- 1:n
  Answer <- mean(X * ((i - 1) / (n - 1)) ^ u)
  return(Answer)
}
#----------------------------
M10v = function(X, v) {
  #v can ONLY be 0 or 1
  n <- length(X)
  X <- sort(X)
  i <- 1:n
  Answer <- mean(X * (1 - ((i - 1) / (n - 1))) ^ v)
  return(Answer)
}
#----------------------------
betahat.PWM <- function(x) {
  M100 <- M10v(x, 0)
  M101 <- M10v(x, 1)
  answer <- (2 * M101 - M100) / (4 * M101 - M100)
  return(answer)
}
#----------------------------
sigmahat.PWM <- function(x) {
  M100 <- M10v(x, 0)
  M101 <- M10v(x, 1)
  answer <- (2 * M100 * M101) / (M100 - 4 * M101)
  return(answer)
}
#----------------------------
PWM <- function(x) {
  n <- length(x)
  x <- sort(x)
  bb <- betahat.PWM(x)
  ss <- sigmahat.PWM(x)
  return(list(beta = bb, sigma = ss))
}
#--------------------------------------
#Airborne Communication Transciever Data
aircomtra<- c(0.2, 0.3, 0.5, 0.5, 0.5, 0.5, 0.6, 0.6 ,0.7 ,0.7 ,0.7 ,0.8,
              0.8, 1.0, 1.0, 1.0, 1.0, 1.1, 1.3, 1.5 ,1.5 ,1.5 ,1.5 ,2.0,
              2.0, 2.2, 2.5, 2.7, 3.0, 3.0, 3.3, 3.3 ,4.0 ,4.0 ,4.5 ,4.7,
              5.0, 5.4, 5.4, 7.0, 7.5, 8.8, 9.0, 10.3, 22.0, 24.5)

graphsnstuff <- function(x,B=1,PvalAIC="p",namex="x",bsCI=FALSE,alpha=0.05,printPDF=TRUE,xlimmax=5){
  tim<-format(Sys.time(), "%Y-%m-%d.%Hh%Mm%Ss")
  varnm<-deparse(substitute(x))
  tim <- paste0(varnm,".",tim)
  
  est <- list()
  
  est[[paste0("MLE",varnm)]]           <- MLEx           <- MLE(x)
  est[[paste0("MLE.b",varnm)]]         <- MLE.bx         <- MLE.b(x)
  est[[paste0("LME",varnm)]]           <- LMEx           <- LME(x)
  est[[paste0("MDE.CvM",varnm)]]       <- MDE.CvMx       <- MDE.CvM(x)
  est[[paste0("MDE.LS",varnm)]]        <- MDE.LSx        <- MDE.LS(x)
  est[[paste0("MDE.Phi.chisq",varnm)]] <- MDE.Phi.chisqx <- MDE.Phi.chisq(x)
  est[[paste0("MDE.Phi.tv",varnm)]]    <- MDE.Phi.tvx    <- MDE.Phi.tv(x)
  est[[paste0("MDE.Phi.kl",varnm)]]    <- MDE.Phi.klx    <- MDE.Phi.kl(x)
  est[[paste0("MME",varnm)]]           <- MMEx           <- MME(x)
  est[[paste0("PWM",varnm)]]           <- PWMx           <- PWM(x)
  
  EDF <- ecdf(x)
  
  
  #LOMAX
  KS.stat <- NULL
  KS.stat[1] <- max(abs(EDF(x) - pLomax(x,MLEx$beta,MLEx$sigma)))
  KS.stat[2] <- max(abs(EDF(x) - pLomax(x,MDE.LSx$beta,MDE.LSx$sigma)))
  
  #WEIBULL
  weib.est   <- eweibull(x)$parameters
  weib.shape <- weib.est["shape"]
  weib.scale <- weib.est["scale"]
  KS.stat[3] <- max(abs(EDF(x) - pweibull(x,shape=weib.shape,scale=weib.scale)))
  
  #GAMMA
  gamm.est   <- MASS::fitdistr(x, "gamma", start=list(shape=1, scale=1))$estimate
  gamm.shape <- gamm.est["shape"]
  gamm.scale <- gamm.est["scale"]
  KS.stat[4] <- max(abs(EDF(x) - pgamma(x,shape=gamm.shape,scale=gamm.scale)))
  
  
  if(bsCI==TRUE){
    pb <- txtProgressBar(min = 0, max = B, initial = 0,style=3) 
    
    b <- 1
    BS.stat <- array(dim=c(B,2,9))
    while(b<=B){
      xstar              <- sample(x,replace=TRUE)
      MLExstar           <- try(MLE(xstar),silent=TRUE)
      MLE.bxstar         <- try(MLE.b(xstar),silent=TRUE)
      LMExstar           <- try(LME(xstar),silent=TRUE)
      MDE.CvMxstar       <- try(MDE.CvM(xstar),silent=TRUE)
      MDE.LSxstar        <- try(MDE.LS(xstar),silent=TRUE)
      MDE.Phi.chisqxstar <- try(MDE.Phi.chisq(xstar),silent=TRUE)
      MDE.Phi.tvxstar    <- try(MDE.Phi.tv(xstar),silent=TRUE)
      MDE.Phi.klxstar    <- try(MDE.Phi.kl(xstar),silent=TRUE)
      MMExstar           <- try(MME(xstar),silent=TRUE)

      if(!any(is.character( MLExstar          ),
              is.character( MLE.bxstar        ),
              is.character( LMExstar          ),
              is.character( MDE.CvMxstar      ),
              is.character( MDE.LSxstar       ),
              is.character( MDE.Phi.chisqxstar),
              is.character( MDE.Phi.tvxstar   ),
              is.character( MDE.Phi.klxstar   ),
              is.character( MMExstar          )
      )){
        
        BS.stat[b,,1]    <- unlist(MLExstar[c("beta","sigma")])          
        BS.stat[b,,2]    <- unlist(MLE.bxstar[c("beta","sigma")])          
        BS.stat[b,,3]    <- unlist(LMExstar[c("beta","sigma")])          
        BS.stat[b,,4]    <- unlist(MDE.CvMxstar[c("beta","sigma")])          
        BS.stat[b,,5]    <- unlist(MDE.LSxstar[c("beta","sigma")])          
        BS.stat[b,,6]    <- unlist(MDE.Phi.chisqxstar[c("beta","sigma")])          
        BS.stat[b,,7]    <- unlist(MDE.Phi.tvxstar[c("beta","sigma")])          
        BS.stat[b,,8]    <- unlist(MDE.Phi.klxstar[c("beta","sigma")])          
        BS.stat[b,,9]    <- unlist(MMExstar[c("beta","sigma")])          

        b <- b+1
        setTxtProgressBar(pb,b)
      }
    }
    close(pb)
    est1 <- data.frame(
           MLE=           unlist(MLEx[c("beta","sigma")]),          
           MLE.b=         unlist(MLE.bx[c("beta","sigma")]),
           LME=           unlist(LMEx[c("beta","sigma")]),
           MDE.CvM=       unlist(MDE.CvMx[c("beta","sigma")]),
           MDE.LS=        unlist(MDE.LSx[c("beta","sigma")]),
           MDE.Phi.chisq= unlist(MDE.Phi.chisqx[c("beta","sigma")]),
           MDE.Phi.tv=    unlist(MDE.Phi.tvx[c("beta","sigma")]),
           MDE.Phi.kl=    unlist(MDE.Phi.klx[c("beta","sigma")]),
           MME=           unlist(MMEx[c("beta","sigma")])
           )
    betaCIs  <- as.data.frame(lapply(data.frame(BS.stat[,1,]),sort))[c(floor((B+1)*(alpha/2)),floor((B+1)*(1-alpha/2))),]
    sigmaCIs <- as.data.frame(lapply(data.frame(BS.stat[,2,]),sort))[c(floor((B+1)*(alpha/2)),floor((B+1)*(1-alpha/2))),]
    
    names(betaCIs) <- names(sigmaCIs) <- names(est1)
    
    est1 <- rbind(est1,betaCIs,sigmaCIs)[c(1,3,4,2,5,6),]
    dimnames(est1) <- list(c("beta",paste0("beta.",(100*(1-alpha)),"%.bsCI:",100*(alpha/2)),paste0("beta.",(100*(1-alpha)),"%.bsCI:",100*(1-alpha/2)),"sigma",paste0("sigma.",(100*(1-alpha)),"%.bsCI:",100*(alpha/2)),paste0("sigma.",(100*(1-alpha)),"%.bsCI:",100*(1-alpha/2))),names(est1))
    
    print(est1)
    
  }  
  
  if(tolower(PvalAIC)=="p"){
    
    b <- 1
    BS.stat <- matrix(nrow=B,ncol=4)
    while(b<=B){
      xstar            <- sample(x,replace=TRUE)
      EDFstar          <- ecdf(xstar)(xstar)
      MLExstar         <- try(MLE(xstar),silent=TRUE)
      MDE.LSxstar      <- try(MDE.LS(xstar),silent=TRUE)
      weib.est.star    <- try(eweibull(xstar)$parameters,silent=TRUE)
      gamm.est.star   <- MASS::fitdistr(xstar, "gamma", start=list(shape=1, scale=1))$estimate
      
      if(!any(is.character(MLExstar),
              is.character(MDE.LSxstar),
              is.character(weib.est.star),
              is.character(gamm.est.star)
      )){
        #LOMAX - MLE
        BS.stat[b,1]    <- max( abs(EDFstar - pLomax(xstar,MLExstar$beta,MLExstar$sigma) ))
        
        #LOMAX - MDE.SD
        BS.stat[b,2]    <- max( abs(EDFstar - pLomax(xstar,MDE.LSxstar$beta,MDE.LSxstar$sigma) ))
        
        #WEIBULL
        weib.shape.star <- weib.est.star["shape"]
        weib.scale.star <- weib.est.star["scale"]
        BS.stat[b,3]    <- max( abs(EDFstar - pweibull(xstar,shape=weib.shape.star,scale=weib.scale.star)))
        
        #GAMMA
        gamm.shape.star <- gamm.est.star["shape"]
        gamm.scale.star <- gamm.est.star["scale"]
        BS.stat[b,4]    <- max( abs(EDFstar - pgamma(xstar,shape=gamm.shape.star,scale=gamm.scale.star)))
        
        b <- b+1
      }
    }
    
    pvalue <- numeric(4)
    pvalue[1] <- mean(BS.stat[,1] > KS.stat[1])
    pvalue[2] <- mean(BS.stat[,2] > KS.stat[2])
    pvalue[3] <- mean(BS.stat[,3] > KS.stat[3])
    pvalue[4] <- mean(BS.stat[,4] > KS.stat[4])
    
    pvalue.char <- as.character(pvalue)
    pvalue.char <- paste0(" (",pvalue.char,")")
    pvalue.char <- c(pvalue.char,"")
    
    my.legend <- paste0(c("Lomax - MLE","Lomax - MDE.SD","Weibull","Gamma","EDF"),pvalue.char)
    my.col    <- c("red","forestgreen","orange2","dodgerblue3","mediumpurple1")
    my.lty    <- c(1,2,4,5,1)
    my.lwd    <- c(1.5,1.5,1.5,1.5,2)
    
    xseq <- seq(min(x)*0.9,max(x)*1.1,length=1000)
    
    if(printPDF)
      pdf(file=paste0("C:/Temp/FittedDists.CDF.",tim,".pdf"),width=8,height=4)
    plot(ecdf(x),main="",xlab=namex,ylab=expression(P(X <= x))             ,lty=my.lty[5],lwd=my.lwd[5],col=my.col[5], xlim=c(0,xlimmax))
    legend(x="bottomright",legend=my.legend                                ,lty=my.lty   ,lwd=my.lwd   ,col=my.col)
    lines(xseq,             pLomax(xseq,    MLEx$beta,     MLEx$sigma)     ,lty=my.lty[1],lwd=my.lwd[1],col=my.col[1])
    lines(xseq,             pLomax(xseq,    MDE.LSx$beta,  MDE.LSx$sigma)  ,lty=my.lty[2],lwd=my.lwd[2],col=my.col[2])
    lines(xseq,           pweibull(xseq,shape=weib.shape,scale=weib.scale) ,lty=my.lty[3],lwd=my.lwd[3],col=my.col[3])
    lines(xseq,             pgamma(xseq,shape=gamm.shape,scale=gamm.scale) ,lty=my.lty[6],lwd=my.lwd[6],col=my.col[6])
    if(printPDF)
      dev.off()
    
    my.legend <- paste0(c("Lomax - MLE","Lomax - MDE.SD","Weibull","Gamma","BC.KDE"),pvalue.char)
    if(printPDF)
      pdf(file=paste0("C:/Temp/FittedDists.PDF.",tim,".pdf"),width=8,height=4)
    dx        <- data.frame(x=xseq,y=dbckden(xseq,x,bw=bw.nrd0(x)))
    ylimmax   <- max(dx$y)*1.2
    hist(x,main="",xlab=namex,freq=FALSE,ylim=c(0,ylimmax), xlim=c(0,xlimmax))
    legend(x="topright",legend=my.legend                                    ,lty=my.lty   ,lwd=my.lwd   ,col=my.col)
    lines(xseq,             dx$y                                            ,lty=my.lty[5],lwd=my.lwd[5],col=my.col[5])
    lines(xseq,             dLomax(xseq,    MLEx$beta,     MLEx$sigma)      ,lty=my.lty[1],lwd=my.lwd[1],col=my.col[1])
    lines(xseq,             dLomax(xseq,    MDE.LSx$beta,     MDE.LSx$sigma)      ,lty=my.lty[2],lwd=my.lwd[2],col=my.col[2])
    lines(xseq,             dweibull(xseq,shape=weib.shape,scale=weib.scale),lty=my.lty[3],lwd=my.lwd[3],col=my.col[3])
    lines(xseq,             dgamma(xseq,shape=gamm.shape,scale=gamm.scale)  ,lty=my.lty[6],lwd=my.lwd[6],col=my.col[6])
    if(printPDF)
      dev.off()
    answr <- list(estimators=est,p.values=pvalue)
  }else if(PvalAIC == "aic"){
    AIC <- numeric(3)
    #AIC Lomax
    AIC[1]<- 4-2*sum(log(dLomax(x,MLEx$beta,MLEx$sigma)))
    #AIC Weibull
    AIC[2]<- 4-2*sum(log(dweibull(x,shape=weib.shape,scale=weib.scale)))
    #AIC Gamma
    AIC[3]<- 4-2*sum(log(dgamma(x,shape=gamm.shape,scale=gamm.scale)))
    
    names(AIC) <- c("Lomax","Weibull","Gamma")
    
    AIC.char    <- as.character(round(AIC[c(1,1,2,3)],1))
    AIC.char <- paste0(" (",AIC.char,")")
    AIC.char[c(2,5)] <- ""
    
    my.legend <- paste0(c("Lomax - MLE","Lomax - MDE.SD","Weibull","Gamma","EDF"),AIC.char)
    my.col    <- c("red","forestgreen","orange2","dodgerblue3","mediumpurple1")
    my.lty    <- c(1,2,4,5,1)
    my.lwd    <- c(1.5,1.5,1.5,1.5,2)
    
    xseq <- seq(min(x)*0.9,max(x)*1.1,length=1000)
    
    if(printPDF)
      pdf(file=paste0("C:/Temp/FittedDists.CDF.",tim,".pdf"),width=8,height=4)
    plot(ecdf(x),main="",xlab=namex,ylab=expression(P(X <= x))    ,lty=my.lty[5],lwd=my.lwd[5],col=my.col[5], xlim=c(0,xlimmax))
    legend(x="bottomright",legend=my.legend                                ,lty=my.lty   ,lwd=my.lwd   ,col=my.col)
    lines(xseq,             pLomax(xseq,    MLEx$beta,     MLEx$sigma)     ,lty=my.lty[1],lwd=my.lwd[1],col=my.col[1])
    lines(xseq,             pLomax(xseq,    MDE.LSx$beta,  MDE.LSx$sigma)  ,lty=my.lty[2],lwd=my.lwd[2],col=my.col[2])
    lines(xseq,           pweibull(xseq,shape=weib.shape,scale=weib.scale) ,lty=my.lty[3],lwd=my.lwd[3],col=my.col[3])
    lines(xseq,       pgamma(xseq,shape=gamm.shape,scale=gamm.scale)       ,lty=my.lty[4],lwd=my.lwd[4],col=my.col[4])
    if(printPDF)
      dev.off()
    
    my.legend <- paste0(c("Lomax - MLE","Lomax - MDE.SD","Weibull","Gamma","BC.KDE"),AIC.char)
    
    if(printPDF)
      pdf(file=paste0("C:/Temp/FittedDists.PDF.",tim,".pdf"),width=8,height=4)
    dx      <- data.frame(x=xseq,y=dbckden(xseq,x,bw=bw.nrd0(x)))
    ylimmax <- max(dx$y)*1.2
    hist(x,main="",xlab=namex,freq=FALSE,ylim=c(0,ylimmax), xlim=c(0,xlimmax))
    legend(x="topright",legend=my.legend                                    ,lty=my.lty   ,lwd=my.lwd   ,col=my.col)
    lines(xseq,             dx$y                                            ,lty=my.lty[5],lwd=my.lwd[5],col=my.col[5])
    lines(xseq,             dLomax(xseq,    MLEx$beta,     MLEx$sigma)      ,lty=my.lty[1],lwd=my.lwd[1],col=my.col[1])
    lines(xseq,             dLomax(xseq,    MDE.LSx$beta,     MDE.LSx$sigma)      ,lty=my.lty[2],lwd=my.lwd[2],col=my.col[2])
    lines(xseq,             dweibull(xseq,shape=weib.shape,scale=weib.scale),lty=my.lty[3],lwd=my.lwd[3],col=my.col[3])
    lines(xseq,               dgamma(xseq,shape=gamm.shape,scale=gamm.scale),lty=my.lty[4],lwd=my.lwd[4],col=my.col[4])
    if(printPDF)
      dev.off()
    
    answr <- list(estimators=est,AIC=AIC)
  }else{
    warning("Specify PvalAIC correctly (must be either 'p' or 'aic').")
    answr <- NULL
  }
  return(answr)
}

graphsnstuff(aircomtra,B=1000,PvalAIC="aic",namex="Repair times, x",bsCI=TRUE,alpha=0.05,xlimmax=15,printPDF=TRUE)