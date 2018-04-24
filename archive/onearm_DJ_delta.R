library(MASS)
library(nlme)
library(reshape2)
library(parallel)
library(invgamma)
library(MCMCpack)
library(Hmisc)
library(gmodels)
library(mvtnorm)

rm(list = ls())


mu_dj <- c(1, 3)
mu_df <- c(1, 3/2)

sd <- c(1.8, 2.2)
corr <- matrix(c(1.0, 0.8,
                 0.8, 1.0), nrow=2)
vv <- diag(sd) %*% corr %*% diag(sd)
N <- 100
alpha <- .05

# S1
psi1 <- c(1, 0, 0)
psi2 <- c(1, 0, 0)
# S2
# psi1 <- c(1, -0.2, 0)
# psi2 <- c(1, -0.4, 0)
# S3
# psi1 <- c(1, -0.4, 0)
# psi2 <- c(1, -0.8, 0)
# S4
# psi1 <- c(1, 0, -0.2)
# psi2 <- c(1, 0, -0.4)
# S5
# psi1 <- c(2, 0, -0.4)
# psi2 <- c(2, 0, -0.8)
# psi1 <- c(1, -1, 0)
# psi2 <- c(1, -1, 0)

px1 <- psi1 %*% c(1,mu_dj)
px2 <- psi2 %*% c(1,mu_dj)
p1 <- 1/(1+exp(-px1))
p2 <- 1/(1+exp(px2))
p1
p2
S <- 5000

true_dj <- mu_dj[2]
true_df <- (1-p1-p2)*mu_dj[2] + (p1+p2)*mu_df[2]


no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl,varlist=c('S', 'N', 'vv', 'alpha', 
                           'true_dj','true_df',
                           'psi1', 'psi2', 
                           'mu_dj','mu_df'))
clusterSetRNGStream(cl, 12345)
clusterEvalQ(cl, {
  library(MASS)
  library(nlme)
  library(reshape2)
  library(invgamma)
  library(MCMCpack)
  library(Hmisc)
  library(gmodels)
  library(mvtnorm)
})


result <- parSapply(cl, 1:S, function(x){
  out <- c()
  y_dj <- mvrnorm(N, mu_dj, vv)
  y_dj <- data.frame(id=1:N, t1=y_dj[,1], t2=y_dj[,2], arm='trt',stringsAsFactors=FALSE)
  px <- data.frame(p1=as.matrix(cbind(1,y_dj[,2:3])) %*% t(t(psi1)),
                   p2=as.matrix(cbind(1,y_dj[,2:3])) %*% t(t(psi2)))
  pi <- 1/(1+exp(-px))
  pi <- data.frame(p0=1-pi$p1-pi$p2, pi)
  
  ind <- apply(pi, 1, function(x){sample(1:3, 1, replace=T, prob=x)})
  while(length(unique(ind))<3|min(table(ind))<3)
    ind <- apply(pi, 1, function(x){sample(1:3, 1, replace=T, prob=x)})
  y_df <- y_dj
  
  ind_df <- which(ind==2 | ind==3)
  y_df$arm[ind_df] <- 'act'
  y_df$t2[ind_df] <- rnorm(length(ind_df), 
                           mean=mu_df[2]+vv[2,1]/vv[1,1]*(y_df$t1[ind_df]-mu_df[1]),
                           sd=sqrt(vv[2,2]-vv[2,1]/vv[1,1]*vv[1,2]))
  
  
  y_mis <- y_df
  ind_mis <- which(ind==3)
  y_mis$arm[ind_mis] <- 'mis'
  y_mis$t2[ind_mis] <- NA
  
  
  r <- N -length(ind_df)
  q <- N - length(ind_mis)
  
  ## De Jure observed
  y_dj_obs <- y_dj[which(ind==1),]
  theta_dj <- mean(y_dj_obs$t2)
  vartotal_dj <- var(y_dj_obs$t2)
  se_dj <- sqrt(vartotal_dj/r)
  low_dj <- theta_dj - qt(1-alpha/2, r-1)*se_dj
  upp_dj <- theta_dj + qt(1-alpha/2, r-1)*se_dj
  cover_dj <- low_dj<true_dj & true_dj<upp_dj
  
  out <- c(out, theta_dj, se_dj, cover_dj)
  
  ## De Jure missing
  # MMRM
  y_mis_long <- melt(y_mis,id.vars = c('id','arm'))
  colnames(y_mis_long) <- c('id','arm','time','val')
  y_mis_long$delta <- as.numeric(y_mis_long$arm=='act' & y_mis_long$time=='t2')
  y_mis_long <- na.omit(y_mis_long)
  
  m1 <- lme(val~time+delta, 
            random=~1+time|id, 
            method='REML', 
            data=y_mis_long,
            control=lmeControl(returnObject=TRUE))
  est <- estimable(m1,
                   cm=cbind(
                     '(Intercept)'=1,
                     'timet2'=1
                   ),
                   conf.int=0.95)
  est <- as.numeric(est)
  theta_dj_mmrm <- est[1]
  coeff <- summary(m1)$coef$fixed
  se_dj_mmrm <- est[2]
  # se_dj_mmrm <- sqrt(est[2]^2+(coeff[2]+coeff[3])^2/N^3*(q*N+r*q-r^2-q^2))
  low_dj_mmrm <- theta_dj_mmrm - qt(1-alpha/2, N-2)*se_dj_mmrm
  upp_dj_mmrm <- theta_dj_mmrm + qt(1-alpha/2, N-2)*se_dj_mmrm
  cover_dj_mmrm <- low_dj_mmrm<true_dj & true_dj<upp_dj_mmrm
  out <- c(out, theta_dj_mmrm, se_dj_mmrm, cover_dj_mmrm)
  
  # MI
  MI <- 50
  z <- matrix(c(1,1,0,1), nrow=2)
  sigma <- z%*%getVarCov(m1)%*%t(z) + diag(2)*m1$sigma^2
  z_trt <- matrix(c(1,0,0,1,1,0), nrow=2,byrow=T)
  sigma_trt <- z_trt%*%m1$varFix%*%t(z_trt)
  mu_trt <- z_trt %*% summary(m1)$coef$fixed
  theta_dj_mi_j <- vartotal_dj_mi_j <- rep(NA, MI)
  for(j in 1:MI){
    mu_trt_tilde <- rmvnorm(1, mu_trt, sigma_trt)
    sigma_tilde <- riwish(N-1, sigma)*(N-1)
    
    y_mis$t2[ind_df] <- rnorm(length(ind_df),
                              mu_trt_tilde[2]+sigma_tilde[1,2]/sigma_tilde[1,1]*(y_mis$t1[ind_df]-mu_trt_tilde[1]),
                              sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
    
    theta_dj_mi_j[j] <- mean(y_mis$t2)
    vartotal_dj_mi_j[j] <- var(y_mis$t2)/(N-1)
  }
  theta_dj_mi <- mean(theta_dj_mi_j)
  W <- mean(vartotal_dj_mi_j)
  B <- var(theta_dj_mi_j)
  se_dj_mi <-  sqrt(W + (1+1/MI)*B)
  df <- (1+MI/(MI+1)*W/B)*(MI-1)
  low_dj_mi <- theta_dj_mi - qt(1-alpha/2, df)*se_dj_mi
  upp_dj_mi <- theta_dj_mi + qt(1-alpha/2, df)*se_dj_mi
  cover_dj_mi <- low_dj_mi<true_dj & true_dj<upp_dj_mi
  
  out <- c(out, theta_dj_mi, se_dj_mi, cover_dj_mi)
})

stopCluster(cl)

round(rowMeans(result),2)
round(sqrt(apply(result,1,var)),2)
