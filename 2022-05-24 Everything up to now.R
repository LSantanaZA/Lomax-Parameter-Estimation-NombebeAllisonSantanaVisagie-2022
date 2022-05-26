rm(list = ls())
cat("\014")

set.seed(123)

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



#raw moments
muk <- function(x, k) {
  mean(x ^ k)
}


#Lomax distribution function
pLomax <- function(q, beta, sigma) {
  1 - (sigma / (q + sigma)) ^ beta
}
#Lomax quantile function
qLomax <- function(x, beta, sigma) {
  sigma * ((1 - x) ^ (-1 / beta) - 1)
}
#Lomax density function
dLomax <- function(x, beta, sigma) {
  (beta * sigma ^ beta) / ((x + sigma) ^ (beta + 1))
}
#Lomax random number generator
rLomax <- function(n, beta, sigma) {
  U <- runif(n)
  qLomax(U, beta, sigma)
}
# #Log-likelihood function (Lomax) #used in optim()
# llik <- function(param,n=100){
#   beta <- param[1]
#   sigma <- param[2]
#   -(n*log(beta)-n*log(sigma)-(1+beta)*sum(log(1+x/sigma)))
# }

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

#Different phi function choices
phi.kl <-
  function(t) {
    -log(t)
  } #Should be t*log(t), but there are numerical complications
phi.chisq <- function(t) {
  (1 - t) ^ 2
}
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
  #b <- 2 *(mu2 - mu1 ^ 2) / (mu2 - 2 * mu1 ^ 2)
  s <- mu1 * (b - 1)
  ans2 <- c(b, s)
  ans3 <- list(beta = ans2[1], sigma = ans2[2])
  return(ans3)
}
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
MMEsigma <- function(x) {
  n <- length(x)
  mu1 <- mean(x)
  mu2 <- mean(x^2)
  #b <- 2*(mu1 - mu1^2)/(mu2 - 2*mu2^2)
  b <- 2*(mu2 - mu1^2)/(mu2 - 2*mu1^2)
  s <- mu1*(b-1)
  ans2 <- s
  return(ans2)
  
}


# Delta method
DeltaMethodFunc <- function(){
  
  EXd <- function(d,b,s){
    s^d*gamma(b-d)*gamma(1 + d)/gamma(b)
  }
  
  bet <- 5
  sig <- 1
  u <- mu1 <- sig/(bet-1)
  v <- mu2 <- 2*sig^2/((bet-1)*(bet-2))
  
  # mu1 <- u <- 1
  # mu2 <- v <- 10
  # bet <- (2*v - 2*u^2)/(v - 2*u^2)
  # sig <- u*(bet - 1)
  
  n <-200
  
  nabla.g1 <- c(v*(2*u^2 + v)/((v-2*u^2)^2), -2*u^3/((v-2*u^2)^2))
  nabla.g2 <- c(4*u*v/((v-2*u^2)^2), -2*u^2/((v-2*u^2)^2))
  
  # varmu1hat <- (EXd(2,bet,sig) - EXd(1,bet,sig)^2)/n
  # varmu2hat <- (EXd(4,bet,sig) - EXd(2,bet,sig)^2)/n
  
  # 24*sig^4/((bet-1)*(bet-2)*(bet-3)*(bet-4))
  # EXd(4,bet,sig)
  # 
  # 2*sig^2/((bet-1)*(bet-2))
  # EXd(2,bet,sig)
  
  COV <- matrix(c(bet*sig^2/(n*(bet-1)^2*(bet-2)),
                  4*sig^3*bet/(n*(bet-1)^2*(bet-2)*(bet-3)),
                  4*sig^3*bet/(n*(bet-1)^2*(bet-2)*(bet-3)),
                  4*bet*(5*bet-11)*sig^4/(n*(bet-1)^2*(bet-2)^2*(bet-3)*(bet-4))),ncol=2)
  COV
  cov.bet <- t(nabla.g2)%*%COV%*%as.matrix(nabla.g2,ncol=1)
  cov.bet
  cov.sig <- t(nabla.g1)%*%COV%*%as.matrix(nabla.g1,ncol=1)
  cov.sig
  # BETA <- 2*(mu2 - mu1^2)/(mu2 - 2*mu1^2)
  # SIGMA <- mu1*(BETA-1)
  ans <- c(cov.bet,cov.sig)
  return(ans)
}





#Bias correction Method of moments estimator
MME.b <- function(x, B = 100) {
  n <- length(x)
  bstar <- numeric(B)
  sstar <- numeric(B)
  P <- matrix(0,ncol=n,nrow=B)
  bb <- 1
  #for(bb in 1:B){
  while(bb<=B){
    j <- sample(1:n,n,replace=TRUE)
    Xstar <- x[j]
    CV.flag <-  (sqrt((n - 1)/n)*sd(Xstar)/mean(Xstar) < 1)
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
#--------------------------------------------------
#L-moments estimator:
#--------------------------------------------------
#Hosking, J. R. M. (1990). L-moments: Analysis and estimation of distributions using linear combinations of order statistics. Journal of the Royal Statistical Society. Series B (Methodological), 52, 105-124.
l2 = function(X) {
  n <- length(X)
  X <- sort(X)
  nchoose2 <- choose(n, 2)
  j <- sequence(1:n)
  i <- rep(1:n, 1:n)
  Answer <- sum(X[i] - X[j]) / (2 * nchoose2)
  return(Answer)
}
#The individual \hat{\beta} function
betahat.L <- function(x) {
  L2 <- l2(x)
  L1 <- mean(x)
  answer <- L2 / (2 * L2 - L1)
  return(answer)
}
#The individual \hat{\sigma} function
sigmahat.L <- function(x) {
  L2 <- l2(x)
  L1 <- mean(x)
  answer <- (L1 ^ 2 - L1 * L2) / (2 * L2 - L1)
  return(answer)
}
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
    # ans1 <- 1/(12*n) + sum((1-(sigma/(x+sigma))^beta-(2*i-1)/(2*n))^2)
    ans1 <-   1/(12*n) + sum((pLomax(x,beta,sigma)-(2*(1:n)-1)/(2*n))^2)
    return(ans1)
  }
  lme <- LME(x)
  b <- lme$beta
  s <- lme$sigma
  # ans2 <- optim(c(b, s), MD.CvM, method = "Nelder-Mead")$par
  ans2 <- optim(c(b, s), MD.CvM, method = "BFGS")
  #ans2 <- optim(c(b, s), MD.CvM, method = "Nelder-Mead")
  ans3 <- list(beta=ans2$par[1],sigma=ans2$par[2],convergence=ans2$convergence)
  return(ans3)
}
#Minimum Distance Estimator (LS)
MDE.LS <- function(x) {
  n <- length(x)
  MD.LS <- function(param) {
    beta <- param[1]
    sigma <- param[2]
    # ans2 <- (1/n)*sum(((1-(sigma/(x + sigma))^beta)-i/(n+1))^2)
    ans2 <- (1 / n) * sum((pLomax(x, beta, sigma) - (1:n) / (n + 1)) ^ 2)
    return(ans2)
  }
  lme <- LME(x)
  b <- lme$beta
  s <- lme$sigma
  # ans2 <- optim(c(b, s), MD.LS, method = "Nelder-Mead")$par
  #ans2 <- optim(c(b, s), MD.LS, method = "Nelder-Mead")
  ans2 <- optim(c(b, s), MD.LS, method = "BFGS")
  ans3 <- list(beta=ans2$par[1],sigma=ans2$par[2],convergence=ans2$convergence)
  return(ans3)
}
#Minimum Distance Estimator (Phi=KL)
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
  # ans2 <- optim(c(b, s), MD.Phi.kl, method = "Nelder-Mead")$par
  ans2 <- optim(c(b, s), MD.Phi.kl, method = "Nelder-Mead")
  ans3 <- list(beta = ans2$par[1], sigma = ans2$par[2], convergence=ans2$convergence)
  return(ans3)
}
#Minimum Distance Estimator (Phi=ChiSq)
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
  # ans2 <- optim(c(b, s), MD.Phi.chisq, method = "Nelder-Mead")$par
  ans2 <- optim(c(b, s), MD.Phi.chisq, method = "Nelder-Mead")
  # ans3 <- list(beta = ans2[1], sigma = ans2[2])
  ans3 <- list(beta = ans2$par[1], sigma = ans2$par[2], convergence=ans2$convergence)
  return(ans3)
}
#Minimum Distance Estimator (Phi=Total Variation)
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
  # ans2 <- optim(c(b, s), MD.Phi.tv, method = "Nelder-Mead")$par
  ans2 <- optim(c(b, s), MD.Phi.tv, method = "Nelder-Mead")
  # ans3 <- list(beta = ans2[1], sigma = ans2[2])
  ans3 <- list(beta = ans2$par[1], sigma = ans2$par[2], convergence=ans2$convergence)
  return(ans3)
}
#----------------------------------------------------


#GRADIENT Log-likelihood function (Lomax)
grad.llik <- function(param, x) {
  n <- length(x)
  beta <- param[1]
  sigma <- param[2]
  ans1 <- n / beta - sum(log(1 + x / sigma))
  ans2 <- -n / sigma + ((1 + beta) / sigma) * sum(x / (sigma + x))
  ans <- c(ans1, ans2)
  return(ans)
}

#HESSIAN Log-likelihood function (Lomax)
hess.llik <- function(param, x) {
  n <- length(x)
  beta <- param[1]
  sigma <- param[2]
  ans11 <- -(n / beta ^ 2)
  #ans22 <- -n*beta/((beta+2)*sigma^2)
  #ans12 <- n/(sigma*(beta+1))
  ans22 <-
    (n / sigma ^ 2) - ((1 + beta) / sigma ^ 2) * sum(x / (sigma + x)) - ((1 +
                                                                            beta) / sigma) * sum(x / (sigma + x) ^ 2)
  ans12 <- (1 / sigma) * sum(x / (sigma + x))
  ans <- matrix(c(ans11, ans12, ans12, ans22), ncol = 2)
  return(ans)
}

#Maximum Likelihood Estimator
MLE <- function(x){
  ans <- flomax(x)
  ans3 <- list(beta = as.numeric(ans[[1]][1]), sigma = as.numeric(ans[[1]][2]), returnCode = as.numeric(!ans$cvg))
  return(ans3)
}



MLE.old <- function(x) {
  n <- length(x)
  lme <- LME(x)
  b <- lme$beta
  s <- lme$sigma
  #
  # A <- matrix(c(1,0,0,1), 2)
  # B <- matrix(c(0,0),ncol=1)
  #
  #ans<-optim(par=c(b,s),method="L-BFGS-B",lower=c(1e-29,1e-29),fn=llik,n=n)$par
  #ans<-summary(maxLik(llik,start=c(b,s),method="BFGS",constraints=list(ineqA=A, ineqB=B)))$estimate[,"Estimate"]
  #ans<-summary(maxLik(llik,start=c(b,s),method="NR"))$estimate[,"Estimate"]
  
  # ans<-
  #   summary(
  #     maxLik(llik,
  #            start=c(b,s),
  #            grad=grad.llik,
  #            hess=hess.llik,
  #            method="NR",
  #            control=list(tol=1e-40,gradtol=1e-13,reltol=1e-20,iterlim=1e5),
  #            x=x
  #     )
  #   )$estimate[,"Estimate"]
  #
  
  # tol0 <- 1e-40#1e-8
  # gradtol0 <- 1e-13#1e-8
  # reltol0 <- 1e-20#1.5e-08
  # iterlim0 <- 1e5#1e3
  ans <-
    summary(
      maxLik(
        llik,
        start = c(b, s),
        grad = grad.llik,
        hess = hess.llik,
        method = "NR",
        control = list(
          # tol = tol0,
          # gradtol = gradtol0,
          # reltol = reltol0,
          # iterlim = iterlim0,
          qac = "marquardt"
        ),
        x = x
      )
    )
  ansrc <- ans$returnCode
  # cat("\ni=", i, "\tansrc=", ansrc)
  # while ((ansrc == 4)) {
  #   iterlim0 <- iterlim0 * 10
  #   
  #   ans <-
  #     summary(
  #       maxLik(
  #         llik,
  #         start = c(b, s),
  #         grad = grad.llik,
  #         hess = hess.llik,
  #         method = "NR",
  #         control = list(
  #           tol = tol0,
  #           gradtol = gradtol0,
  #           reltol = reltol0,
  #           iterlim = iterlim0,
  #           qac = "marquardt"
  #         ),
  #         x = x
  #       )
  #     )
  #   ansrc <- ans$returnCode
  # }
  ans <- ans$estimate[, "Estimate"]
  #ans<-summary(maxNR(fn=llik,grad=grad.llik,hess=hess.llik,start=c(b,s)))$estimate[,"estimate"]
  # ans<-optim(par=c(b,s),method="L-BFGS-B",lower=c(1e-29,1e-29),fn=llik,n=n)$par
  ans3 <- list(beta = ans[1], sigma = ans[2], returnCode = ansrc)
  return(ans3)
}



# 
# fd <-
#   fitdist(
#     x,
#     distr = "Lomax",
#     method = "mle",
#     start = list(beta = LME(x)$beta, sigma = LME(x)$sigma)
#   )
# me <- MLE(x)
# rbind(unlist(me), fd$estimate)
X <- c(0.08, 2.09, 3.48, 4.87, 6.94, 8.66, 13.11, 23.63, 0.20, 2.23, 3.52, 4.98, 6.97, 9.02, 13.29, 0.40,
       2.26, 3.57, 5.06, 7.09, 9.22, 13.80, 25.74, 0.50, 2.46, 3.64, 5.09, 7.26, 9.47, 14.24, 25.82, 0.51,
       2.54, 3.70, 5.17, 7.28, 9.74, 14.76, 26.31, 0.81, 2.62, 3.82, 5.32, 7.32, 10.06, 14.77, 32.15,
       2.64, 3.88, 5.32, 7.39, 10.34, 14.83, 34.26, 0.90, 2.69, 4.18, 5.34, 7.59, 10.66, 15.96, 36.66,
       1.05, 2.69, 4.23, 5.41, 7.62, 10.75, 16.62, 43.01, 1.19, 2.75, 4.26, 5.41, 7.63, 17.12, 46.12,
       1.26, 2.83, 4.33, 5.49, 7.66, 11.25, 17.14, 79.05, 1.35, 2.87, 5.62, 7.87, 11.64, 17.36, 1.40,
       3.02, 4.34, 5.71, 7.93, 11.79, 18.10, 1.46, 4.40, 5.85, 8.26, 11.98, 19.13, 1.76, 3.25, 4.50,
       6.25, 8.37, 12.02, 2.02, 3.31, 4.51, 6.54, 8.53, 12.03, 20.28, 2.02, 3.36, 6.76, 12.07, 21.73,
       2.07, 3.36, 6.93, 8.65, 12.63, 22.69)

# 
# 
# x <- c(
#   6766,
#   7123,
#   10562,
#   14474,
#   15351,
#   16983,
#   18383,
#   19030,
#   25304,
#   29112,
#   30146,
#   33727,
#   40596,
#   41409,
#   47905,
#   49397,
#   52600,
#   59917,
#   63123,
#   77809,
#   102942,
#   103217,
#   123680,
#   140136,
#   192013,
#   198446,
#   227338,
#   329511,
#   361200,
#   421680,
#   513586,
#   545778,
#   750389,
#   863881,
#   1638000
# )
# x <- x - 5000
# summary(
#   maxLik(
#     llik,
#     start = c(LME(x)$beta, LME(x)$sigma),
#     grad = grad.llik,
#     hess = hess.llik,
#     method = "NR",
#     control = list(
#       tol = tol0,
#       gradtol = gradtol0,
#       reltol = reltol0,
#       iterlim = iterlim0
#     ),
#     x = x
#   )
# )



#
#
# #Maximum Likelihood Estimator
# MLE <- function(x){
#   n <- length(x)
#   lme <- LME(x)
#   b <- lme$beta
#   s <- lme$sigma
#   #
#   # A <- matrix(c(1,0,0,1), 2)
#   # B <- matrix(c(0,0),ncol=1)
#   #
#   #ans<-optim(par=c(b,s),method="L-BFGS-B",lower=c(1e-29,1e-29),fn=llik,n=n)$par
#   #ans<-summary(maxLik(llik,start=c(b,s),method="BFGS",constraints=list(ineqA=A, ineqB=B)))$estimate[,"Estimate"]
#   #ans<-summary(maxLik(llik,start=c(b,s),method="NR"))$estimate[,"Estimate"]
#   ans<-summary(maxLik(llik,grad=grad.llik,hess=hess.llik,start=c(b,s)))$estimate[,"Estimate"]
#   #ans<-summary(maxNR(fn=llik,grad=grad.llik,hess=hess.llik,start=c(b,s)))$estimate[,"estimate"]
#   # ans<-optim(par=c(b,s),method="L-BFGS-B",lower=c(1e-29,1e-29),fn=llik,n=n)$par
#   ans3 <- list(beta=ans[1],sigma=ans[2])
#   return(ans3)
# }
#
# #-------------------------
# x <- rLomax(100,2,3)
# MLE(x)
# LME(x)
#
# #-------------------------

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
  # ans <- as.vector(solve(K)%*%A%*%as.vector(solve(K)))
  return(list(beta = ans[2], sigma = ans[1]))
}
#---------------------------------------------------------
M1u0 = function(X, u) {
  #u can ONLY be 0 or 1
  n <- length(X)
  X <- sort(X)
  i <- 1:n
  Answer <- mean(X * ((i - 1) / (n - 1)) ^ u)
  return(Answer)
}
M10v = function(X, v) {
  #v can ONLY be 0 or 1
  n <- length(X)
  X <- sort(X)
  i <- 1:n
  Answer <- mean(X * (1 - ((i - 1) / (n - 1))) ^ v)
  return(Answer)
}
betahat.PWM <- function(x) {
  M100 <- M10v(x, 0)
  M101 <- M10v(x, 1)
  answer <- (2 * M101 - M100) / (4 * M101 - M100)
  return(answer)
}
sigmahat.PWM <- function(x) {
  M100 <- M10v(x, 0)
  M101 <- M10v(x, 1)
  answer <- (2 * M100 * M101) / (M100 - 4 * M101)
  return(answer)
}
PWM <- function(x) {
  n <- length(x)
  bb <- betahat.PWM(x)
  ss <- sigmahat.PWM(x)
  return(list(beta = bb, sigma = ss))
}
#--------------------------------------

# n<-20000
# beta <- 8 #8
# sigma <- 15 #3
# orig<-c(beta,sigma)
# names(orig)<-c("beta","sigma")
# x<-rLomax(n,beta,sigma)
# x<-sort(x)
# n<-length(x)
# flomax(x)$estimate
# unlist(MLE(x))
# unlist(MLE.b(x))
# unlist(LME(x))
# unlist(MDE.Phi.chisq(x))
# unlist(MDE.Phi.tv(x))
# unlist(MDE.Phi.kl(x))
# unlist(MME(x))
# unlist(MME.b(x))
# unlist(PWM(x))
# unlist(orig)



starttime <- format(Sys.time(), "%Y%m%d-%Hh%Mm%Ss")

betabeta <-c(1.1,1.5,2,2.1) #c(1.5,2,3,5)
sigsig <- c(2) # c(2) # c(3)
nn <- c(500)
MC <- 10000
B <- 1000
estnames <- c(
  "MLE",
  "MLE.b",
  "LME",
  "MDE.CvM",
  "MDE.LS",
  "MDE.Phi.X2",
  "MDE.Phi.tv",
  "MDE.Phi.kl",
  "MME",
  "MME.b",
  "PWM"
)
parmnames <- c("beta", "sigma")

var.lomax <-
  mu2.lomax <-
  bias.lomax <-
  mse.lomax <-
  mean.lomax <-
  array(
    numeric(
      length(estnames) * length(parmnames) * length(betabeta) * length(sigsig) *
        length(nn)
    ),
    c(
      length(estnames),
      length(parmnames),
      length(betabeta),
      length(sigsig),
      length(nn)
    )
  )

dimnames(mean.lomax) <-
  dimnames(var.lomax) <-
  dimnames(mu2.lomax) <-
  dimnames(bias.lomax) <-
  dimnames(mse.lomax) <-
  list(
    estnames,
    parmnames,
    paste0("beta=", betabeta),
    paste0("sigma=", sigsig),
    paste0("n=", nn)
  )



accuracy.lomax <-
  accuracy2.lomax <-
  array(numeric(
    length(estnames) * length(betabeta) * length(sigsig) * length(nn)
  ),
  c(
    length(estnames),
    1,
    length(betabeta),
    length(sigsig),
    length(nn)
  ))

dimnames(accuracy.lomax) <-
  dimnames(accuracy2.lomax) <-
  list(
    estnames,
    "Measure",
    paste0("beta=", betabeta),
    paste0("sigma=", sigsig),
    paste0("n=", nn)
  )

#FOR THE PROGRESS BAR
i <- 0
pbTOTAL <- length(betabeta) * length(sigsig) * length(nn) * MC
pb <-
  winProgressBar(
    title = "Monte Carlo Progress",
    label = "Calculating...",
    min = 0,
    max = 100,
    initial = 0,
    width = 300
  )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


for (n in nn) {
  for (ss in sigsig) {
    for (bb in betabeta) {
      
      tru.lomax <- matrix(c(rep(bb, length(estnames)), rep(ss, length(estnames))), ncol = 2)
      # tru.lomax <- matrix(c(rep(bb, 9), rep(ss, 9)), ncol = 2)
      dimnames(tru.lomax) <- list(estnames, parmnames)
      
      est.lomax <-
        array(numeric(nrow(tru.lomax) * ncol(tru.lomax) * MC), c(dim(tru.lomax), MC))
      dimnames(est.lomax) <- list(estnames, parmnames, 1:MC)
      
      FF.lomax <-
        array(numeric(nrow(tru.lomax) * 1 * MC), c(nrow(tru.lomax), 1, MC))
      dimnames(FF.lomax) <-
        list(estnames, c("mean.sq.diff.Fthetahat.Ftheta"), 1:MC)
      
      mc <- 1
      while (mc <= MC) {
        x <- rLomax(n, bb, ss)
        x <- sort(x)
        
        CV.flag <-  (sqrt((n - 1)/n)*sd(x)/mean(x) < 1)
        
        # est.lomax[ , ,mc] <- rbind(unlist(MLE(x)),
        #                            unlist(MLE.b(x)),
        #                            unlist(LME(x)),
        #                            unlist(MDE.Phi.chisq(x)),
        #                            unlist(MDE.Phi.tv(x)),
        #                            unlist(MDE.Phi.kl(x)),
        #                            unlist(MME(x)),
        #                            unlist(MME.b(x,B)),
        #                            unlist(PWM(x)))
        
        
        MLEx <- try(unlist(MLE(x)), silent = TRUE)
        # cat("\nMLEx=",MLEx)
        MLE.bx <- try(unlist(MLE.b(x)), silent = TRUE)
        #MLE.bx<-0 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!!!!!!!
        # cat("\nMLE.bx=",MLE.bx)
        LMEx <-
          try(unlist(LME(x)), silent = TRUE)
        # cat("\nLMEx=",LMEx)
        
        MDE.CvMx <-
          try(unlist(MDE.CvM(x)), silent = TRUE)
        # cat("\nLMEx=",LMEx)
        MDE.LSx <-
          try(unlist(MDE.LS(x)), silent = TRUE)
        # cat("\nLMEx=",LMEx)
        
        
        MDE.Phi.chisqx <-
          try(unlist(MDE.Phi.chisq(x)), silent = TRUE)
        # cat("\nMDE.Phi.chisqx=",MDE.Phi.chisqx)
        MDE.Phi.tvx <-
          try(unlist(MDE.Phi.tv(x)), silent = TRUE)
        # cat("\nMDE.Phi.tvx=",MDE.Phi.tvx)
        MDE.Phi.klx <-
          try(unlist(MDE.Phi.kl(x)), silent = TRUE)
        # cat("\nMDE.Phi.klx=",MDE.Phi.klx)
        MMEx <- try(unlist(MME(x)), silent = TRUE)
        # cat("\nMMEx=",MMEx)
        MME.bxB <-
          try(unlist(MME.b(x, B)), silent = TRUE)
        # cat("\nMME.bxB=",MME.bxB)
        PWMx <-
          try(unlist(PWM(x)), silent = TRUE)
        # cat("\nPWMx=",PWMx)
        
        # cat("mc=",mc," T/F=",
        #   !is.character(MLEx),
        #   !is.character(MLE.bx),
        #   !is.character(LMEx),
        #   !is.character(MDE.CvMx),
        #   !is.character(MDE.LSx),
        #   !is.character(MDE.Phi.chisqx),
        #   !is.character(MDE.Phi.tvx),
        #   !is.character(MDE.Phi.klx),
        #   !is.character(MMEx),
        #   !is.character(MME.bxB),
        #   !is.character(PWMx),"\n")
        
        
        # cat("CV.flag=\t\t",CV.flag,"\n",
        #   "MLEx=\t\t\t",(MLEx),"\n",
        #     "MLE.bx=\t\t",(MLE.bx),"\n",
        #     "LMEx=\t\t\t",(LMEx),"\n",
        #     "MDE.CvMx=\t\t",(MDE.CvMx),"\n",
        #     "MDE.LSx=\t\t",(MDE.LSx),"\n",
        #     "MDE.Phi.chisqx=\t",(MDE.Phi.chisqx),"\n",
        #     "MDE.Phi.tvx=\t\t",(MDE.Phi.tvx),"\n",
        #     "MDE.Phi.klx=\t\t",(MDE.Phi.klx),"\n",
        #     "MMEx=\t\t\t",(MMEx),"\n",
        #     "MME.bxB=\t\t",(MME.bxB),"\n",
        #     "PWMx=\t\t\t",(PWMx),"\n")
        
        # flomax(x,scaleData=FALSE)$estimate
        # flomax(x,scaleData=FALSE)$cvg
        
        #if(CV.flag){
        #if(MMEx["beta"]<0 | MMEx["sigma"]<0){
        # if(MLEx[["returnCode"]]%in% c(3,4,5,6,7,9,100)){
        #   x=x
        # }
        
        
        
        
        if (!any(
          is.character(MLEx),
          is.character(MLE.bx),
          is.character(LMEx),
          is.character(MDE.CvMx),
          is.character(MDE.LSx),
          is.character(MDE.Phi.chisqx),
          is.character(MDE.Phi.tvx),
          is.character(MDE.Phi.klx),
          is.character(MMEx),
          is.character(MME.bxB),
          is.character(PWMx))){
          if(!any(CV.flag,
                  #MME.bxB["sigma"] < 0,
                  LMEx["sigma"] < 0,
                  LMEx["beta"] < 0,
                  MLEx[["returnCode"]] %in% c(3,4,5,6,7,9,100),
                  MDE.CvMx[["convergence"]] >0,
                  MDE.LSx[["convergence"]] >0,
                  MDE.Phi.chisqx[["convergence"]] >0, 
                  MDE.Phi.klx[["convergence"]] >0,
                  MDE.Phi.tvx[["convergence"]] >0
          )) {
            
            #----------------------------------
            # if(MDE.CvMx[1] > 1e3 | MDE.CvMx[2] > 1e3){
            #   banana<-1
            # }
            #----------------------------------  
            
            est.lomax[, , mc] <- rbind(
              MLEx[1:2],
              MLE.bx,
              LMEx,
              MDE.CvMx[1:2],
              MDE.LSx[1:2],
              MDE.Phi.chisqx[1:2],
              MDE.Phi.tvx[1:2],
              MDE.Phi.klx[1:2],
              MMEx,
              MME.bxB,
              PWMx
            )
            
            
            
            FF.lomax[, 1 , mc] <-
              rbind(
                mean((
                  pLomax(x, MLEx["beta"], MLEx["sigma"]) - pLomax(x, bb, ss)
                ) ^ 2),
                mean((
                  pLomax(x, MLE.bx["beta"], MLE.bx["sigma"]) - pLomax(x, bb, ss)
                ) ^ 2),
                mean((
                  pLomax(x, LMEx["beta"], LMEx["sigma"]) - pLomax(x, bb, ss)
                ) ^ 2),
                #
                mean((
                  pLomax(x, MDE.CvMx["beta"], MDE.CvMx["sigma"]) - pLomax(x, bb, ss)
                ) ^ 2),
                mean((
                  pLomax(x, MDE.LSx["beta"], MDE.LSx["sigma"]) - pLomax(x, bb, ss)
                ) ^ 2),
                #
                mean((
                  pLomax(x, MDE.Phi.chisqx["beta"], MDE.Phi.chisqx["sigma"]) - pLomax(x, bb, ss)
                ) ^ 2),
                mean((
                  pLomax(x, MDE.Phi.tvx["beta"], MDE.Phi.tvx["sigma"]) - pLomax(x, bb, ss)
                ) ^ 2),
                mean((
                  pLomax(x, MDE.Phi.klx["beta"], MDE.Phi.klx["sigma"]) - pLomax(x, bb, ss)
                ) ^ 2),
                mean((
                  pLomax(x, MMEx["beta"], MMEx["sigma"]) - pLomax(x, bb, ss)
                ) ^ 2),
                mean((
                  pLomax(x, MME.bxB["beta"], MME.bxB["sigma"]) - pLomax(x, bb, ss)
                ) ^ 2),
                mean((
                  pLomax(x, PWMx["beta"], PWMx["sigma"]) - pLomax(x, bb, ss)
                ) ^ 2)
              )
            
            # print("pLOMAX:------------vvvv")
            # print(pLomax(x,MME.bxB["beta"],MME.bxB["sigma"]))
            # print("BETA:------------vvvv")
            # print(MME.bxB["beta"])
            # print("SIGMA:------------vvvv")
            # print(MME.bxB["sigma"])
            # print("X VALUES:------------vvvv")
            # print(x)
            
            #UPDATE THE PROGRESS BAR
            i <- i + 1 / pbTOTAL
            info <- sprintf("%g%% done", round(i, 3) * 100)
            setWinProgressBar(pb, i * 100, label = info)
            
            mc <- mc + 1
            
          } #end if checking for valid convergence
        }#end if checking for non-character values
        
        
      } #END MC LOOP
      
      tmpmean <- apply(est.lomax, c(1, 2), mean)
      mean.lomax[, , which(bb == betabeta), which(ss == sigsig), which(n ==nn)] <- tmpmean
      
      tmpvar <- apply(est.lomax, c(1, 2), var)
      var.lomax[, , which(bb == betabeta), which(ss == sigsig), which(n ==nn)]  <- tmpvar
      
      tmpmu2 <- apply(est.lomax, c(1, 2), muk, k = 2)
      mu2.lomax[, , which(bb == betabeta), which(ss == sigsig), which(n ==nn)]  <- tmpmu2
      
      bias.lomax[, , which(bb == betabeta), which(ss == sigsig), which(n ==nn)] <- tmpmean - tru.lomax
      
      mse.lomax[, , which(bb == betabeta), which(ss == sigsig), which(n ==nn)]  <- tmpmu2 - 2 * tru.lomax * tmpmean + tru.lomax ^ 2
      
      big.tru.lomax <- array(tru.lomax , dim(est.lomax))
      
      tmp1 <- (est.lomax - big.tru.lomax) ^ 2
      tmp2 <- apply(tmp1, c(1, 3), sum)
      
      accuracy.lomax[, , which(bb == betabeta), which(ss == sigsig), which(n ==nn)]  <-t(apply(tmp2, 1, mean))
      
      tmpFFmean <- apply(FF.lomax, c(1, 2), mean)
      accuracy2.lomax[, , which(bb == betabeta), which(ss == sigsig), which(n ==nn)] <-tmpFFmean
      
      # accuracy.mean <- apply(accuracy.lomax,c(1,2),mean)
    }# END BETABETA
  }#END SIGSIG
}#END NN

#CLOSE THE PROGRESS BAR
close(pb)


mean.lomax
var.lomax
bias.lomax
mse.lomax
accuracy.lomax
accuracy2.lomax



if (!dir.exists("C:/Temp"))
  dir.create("C:/Temp") #Creates the folder if it doesn't already exist.

print(xtable(data.frame(ftable(mean.lomax)), type = "latex"), file = 
        paste0("C:/Temp/",starttime,"-Output.mean.LaTeX.tex"))
print(xtable(data.frame(ftable(var.lomax)), type = "latex"), file = 
        paste0("C:/Temp/",starttime,"-Output.var.LaTeX.tex"))
print(xtable(data.frame(ftable(bias.lomax)), type = "latex"), file = 
        paste0("C:/Temp/",starttime,"-Output.bias.LaTeX.tex"))
print(xtable(data.frame(ftable(mse.lomax)), type = "latex"), file = 
        paste0("C:/Temp/",starttime,"-Output.mse.LaTeX.tex"))
print(xtable(data.frame(ftable(accuracy.lomax)), type = "latex"), file =
        paste0("C:/Temp/",starttime,"-Output.accuracy.LaTeX.tex"))
print(xtable(data.frame(ftable(accuracy2.lomax)), type = "latex"), file =
        paste0("C:/Temp/",starttime,"-Output.accuracy2.LaTeX.tex"))
write.csv(data.frame(ftable(mean.lomax)), file = 
            paste0("C:/Temp/",starttime,"-Output.mean.CSV.csv"))
write.csv(data.frame(ftable(var.lomax)), file = 
            paste0("C:/Temp/",starttime,"-Output.var.CSV.csv"))
write.csv(data.frame(ftable(bias.lomax)), file = 
            paste0("C:/Temp/",starttime,"-Output.bias.CSV.csv"))
write.csv(data.frame(ftable(mse.lomax)), file = 
            paste0("C:/Temp/",starttime,"-Output.mse.CSV.csv"))
write.csv(data.frame(ftable(accuracy.lomax)), file = 
            paste0("C:/Temp/",starttime,"-Output.accuracy.CSV.csv"))
write.csv(data.frame(ftable(accuracy2.lomax)), file = 
            paste0("C:/Temp/",starttime,"-Output.accuracy2.CSV.csv"))


dput(est.lomax,file=
       paste0("C:/Temp/",starttime,"-est.lomax.DUMP.txt"))


write.csv(est.lomax,file=
            paste0("C:/Temp/",starttime,"-est.lomax.EVERYTHING.csv"))
est.lomax0 <- t(read.csv(file=
                           paste0("C:/Temp/",starttime,"-est.lomax.EVERYTHING.csv"),header=TRUE))
write.csv(est.lomax0,file=
            paste0("C:/Temp/",starttime,"-est.lomax.EVERYTHING.csv"))


sink(file=
       paste0("C:/Temp/",starttime,"-est.lomax.PRINT.txt"))
print(est.lomax)
sink()



# pLomax(x,MME.bxB["beta"],MME.bxB["sigma"])
# accuracy2.lomax["MME.b",,,,]
tbl <- ""
for (n in nn) {
  for (ss in sigsig) {
    for (bb in betabeta) {
      tbl <- paste0(
        "% START DATE AND TIME=",
        starttime,
        "
      \\begin{table}[htbp]
  \\centering
  \\caption{$\\beta=",
        bb,
        ",  \\sigma=",
        ss,
        "$ $n=",
        n,
        "$}
  \\resizebox{\\columnwidth}{!}{
    \\begin{tabular}{|cc|r|r|r|r|r|r|r|r|r|r|r|}
      %\\cline{2-13}    \\multicolumn{1}{r|}{} &       & \\multicolumn{11}{c|}{$\\beta=",
        bb,
        ",  \\sigma=",
        ss,
        "$ $n=",
        n,
        "$} \\\\
\\cline{2-13}    \\multicolumn{1}{r|}{} & \\textbf{METHODS} & \\multicolumn{1}{l|}{\\textbf{MLE}} & \\multicolumn{1}{l|}{\\textbf{MLE.b}} & \\multicolumn{1}{l|}{\\textbf{LME}} & \\multicolumn{1}{l|}{\\textbf{MDE.CvM}} & \\multicolumn{1}{l|}{\\textbf{MDE.LS}} & \\multicolumn{1}{l|}{\\textbf{MDE.Phi.X2}} & \\multicolumn{1}{l|}{\\textbf{MDE.Phi.tv}} &\\multicolumn{1}{l|}{\\textbf{MDE.Phi.kl}} & \\multicolumn{1}{l|}{\\textbf{MME}} & \\multicolumn{1}{l|}{\\textbf{MME.b}} & \\multicolumn{1}{l|}{\\textbf{PWM}} \\\\
    \\hline
\\multirow{2}[4]{*}{Mean} & \\multicolumn{1}{|c|}{$\\widehat\\beta$} & ",
        paste0(round(mean.lomax[, 1, which(bb == betabeta), which(ss ==
                                                                    sigsig), which(n == nn)], 4), collapse = "&"),
        "
\\\\
\\cline{2-11}          & \\multicolumn{1}{|c|}{$\\widehat\\sigma$} &
",
        paste0(round(mean.lomax[, 2, which(bb == betabeta), which(ss ==
                                                                    sigsig), which(n == nn)], 4), collapse = "&"),
        "
\\\\
\\hline
\\multirow{2}[4]{*}{Var} & \\multicolumn{1}{|c|}{$\\widehat\\beta$} & ",
        paste0(round(var.lomax[, 1, which(bb == betabeta), which(ss ==
                                                                   sigsig), which(n == nn)], 4), collapse = "&"),
        "
\\\\
\\cline{2-11}          & \\multicolumn{1}{|c|}{$\\widehat\\sigma$} &
",
        paste0(round(var.lomax[, 2, which(bb == betabeta), which(ss ==
                                                                   sigsig), which(n == nn)], 4), collapse = "&"),
        "
\\\\
\\hline
\\multirow{2}[4]{*}{Bias} & \\multicolumn{1}{|c|}{$\\widehat\\beta$} & ",
        paste0(round(bias.lomax[, 1, which(bb == betabeta), which(ss ==
                                                                    sigsig), which(n == nn)], 4), collapse = "&"),
        "
\\\\
\\cline{2-11}          & \\multicolumn{1}{|c|}{$\\widehat\\sigma$} &
",
        paste0(round(bias.lomax[, 2, which(bb == betabeta), which(ss ==
                                                                    sigsig), which(n == nn)], 4), collapse = "&"),
        "
\\\\
\\hline
\\multirow{2}[4]{*}{MSE} & \\multicolumn{1}{|c|}{$\\widehat\\beta$} & ",
        paste0(round(mse.lomax[, 1, which(bb == betabeta), which(ss ==
                                                                   sigsig), which(n == nn)], 4), collapse = "&"),
        "
\\\\
\\cline{2-11}          & \\multicolumn{1}{|c|}{$\\widehat\\sigma$} &
",
        paste0(round(mse.lomax[, 2, which(bb == betabeta), which(ss ==
                                                                   sigsig), which(n == nn)], 4), collapse = "&"),
        "
\\\\
\\hline
\\multicolumn{2}{|c|}{Accuracy (Euclid Dist.)} & ",
        paste0(round(accuracy.lomax[, 1, which(bb == betabeta), which(ss ==
                                                                        sigsig), which(n == nn)], 4), collapse = "&")
        ,
        "\\\\
\\hline
\\multicolumn{2}{|c|}{Accuracy (DF sq.diff.)} & ",
        paste0(round(accuracy2.lomax[, 1, which(bb == betabeta), which(ss ==
                                                                         sigsig), which(n == nn)], 4), collapse = "&")
        ,
        "\\\\
\\hline
\\end{tabular}
}
\\end{table}"
      )
      cat(
        tbl,
        file = paste0(
          "C:/Temp/",
          starttime,
          "-Full.beta=",
          bb,
          ".sig=",
          ss,
          ".n=",
          n,
          ".LaTeX.tex"
        )
      )
    }
  }
}
