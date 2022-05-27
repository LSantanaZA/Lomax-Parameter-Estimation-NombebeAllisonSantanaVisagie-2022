rm(list = ls())
cat("\014")

set.seed(123)

#----------------------------
#Install and load packages:
#----------------------------
#install.packages("Renext")
#library("Renext")
#install.packages("maxLik")
#library("maxLik")
#install.packages("xtable")
#library(xtable)
#install.packages("miscTools")
#library("miscTools")
#install.packages("plm", dependencies=TRUE)
#library("plm")
instll <- c("xtable", "maxLik", "Renext","miscTools","plm")
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
#Bias adjusted Method of moments estimator
#MME.b.old <- function(x, B = 1000) {
#n <- length(x)
#mme <- MME(x)
#datX <- matrix(sample(x, n * B, replace = TRUE), ncol = n)
#b <- 2 * mme$beta - mean(apply(datX, 1, MMEbeta))
#s <- 2 * mme$sigma - mean(apply(datX, 1, MMEsigma))
#ans3 <- list(beta = b, sigma = s)
#return(ans3)
# }
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
  bb <- betahat.PWM(x)
  ss <- sigmahat.PWM(x)
  return(list(beta = bb, sigma = ss))
}
#--------------------------------------


x <- c(
1.58 , 1.65 , 1.73 , 1.81 , 1.88 , 1.96 , 2.04 , 2.12 ,
2.19 , 2.27 , 2.35 , 2.42 , 2.70 , 2.90 , 3.10 , 3.30 ,
3.75 , 4.00 , 4.25 , 4.70 , 4.90 , 5.10 , 5.30 , 5.70 ,
5.90 , 6.10 , 6.30 , 7.83 , 8.17 , 9.00 , 15.00 , 17.00 ,
22.00 , 23.00 , 23.83 , 24.17 , 25.00 , 27.00 , 32.00 , 43.00)
(MLEx <- MLE(x))
(MLE.bx <- MLE.b(x))
(LMEx <- LME(x))
(MDE.CvMx <- MDE.CvM(x))
(MDE.LSx <- MDE.LS(x))
(MDE.Phi.chisqx <- MDE.Phi.chisq(x))
(MDE.Phi.tvx <- MDE.Phi.tv(x))
(MDE.Phi.klx <- MDE.Phi.kl(x))
(MMEx <- MME(x))
(PWMx <- PWM(x))