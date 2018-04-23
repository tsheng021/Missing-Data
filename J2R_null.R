
library(parallel)

rm(list = ls())


mu_cc <- c(10, 20)
mu_ci <- c(10, 20)
mu_ms <- c(10, 20)

vv <- matrix(c(1, .5, .5, 1),nrow=2)
N <- 500
pi <- .5
pi_c <- .3

S <- 10000

# Hypothetical complete
mu2hat <- var2hat <- rep(NA, S)
# No imputation
mu2hat_no <- var2hat_no <- rep(NA, S)
# MMRM
mu2hat_mmrm <- var2hat_mmrm <- rep(NA, S)


no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl,varlist=c('S','pi_c','pi','N', 'vv', 'mu_cc', 'mu_ci', 'mu_ms'))
clusterSetRNGStream(cl, 12345)
clusterEvalQ(cl, {
  library(MASS)
  library(nlme)
  library(reshape2)
  library(parallel)
})


result <- parSapply(cl, 1:S, function(x){
  ## Estimand J2Z ##
  n <- sample(1:3,N,replace=TRUE,prob=c(pi,pi_c,1-pi-pi_c))  # partition size
  n_cc <- sum(n==1)
  n_ci <- sum(n==2)
  n_ms <- sum(n==3)
  
  ## Complete Case ##
  if (n_cc>0){
    y_cc <- data.frame(id=as.character(1:n_cc), mvrnorm(n_cc, mu_cc, vv), group = 'complete')
  }
  colnames(y_cc) <- c('id', 't1', 't2', 'group')
  
  ## Collected Intercurrent ##
  if (n_ci>0){
    y_ci <- data.frame(id=as.character((n_cc+1):(n_cc+n_ci)), mvrnorm(n_ci, mu_ci, vv), group='switch')
  }
  colnames(y_ci) <- c('id', 't1', 't2', 'group')
  
  ## Missing Case ##
  if(n_ms>0){
    y_ms <- data.frame(id=as.character((n_cc+n_ci+1):N) ,mvrnorm(n_ms, mu_ms, vv), group='miss')
  }
  colnames(y_ms) <- c('id', 't1', 't2', 'group')
  
  # The hypothetical complete data #
  y_hyp <- rbind(y_cc, y_ci, y_ms)
  # The observed complete data #
  y_obs <- subset(y_hyp, group!='miss')
  
  mu2hat <- mean(y_hyp$t2)
  var2hat <- var(y_hyp$t2)/N
  
  ## No Imputation ##
  r <- n_cc
  q <- n_cc + n_ci
  mu2hat_no <- r/N*mean(y_ci$t2) + (N-r)/N*mean(y_ci$t2)
  var2hat_no <- var(y_obs$t2)/N
  
  ## MMRM ##
  y_obs_long <- melt(y_obs)
  colnames(y_obs_long) <- c('id','group','time','y')
  m <- lme(y~time*group, 
           random=~1+time*group|id, 
           method='REML', 
           data=y_obs_long,
           control=lmeControl(returnObject=TRUE))
  mu_switch <- sum(coef(summary(m))[,1])
  mu_complete <- sum(coef(summary(m))[1:2,1])
  mu2hat_mmrm <- r/N*mu_complete+(N-r)/N*mu_switch
  var2hat_mmrm <- (m$sigma^2+sum(getVarCov(m)))/N
  
  ## MI Rubin ##
  
  
  
  
  
  r <- c(mu2hat, var2hat, 
         mu2hat_no, var2hat_no, 
         mu2hat_mmrm, var2hat_mmrm)
})

stopCluster(cl)

rowMeans(result)
apply(result,1,var) 

# 
# mean(mu2hat)
# var(mu2hat)
# mean(var2hat)
# 
# mean(mu2hat_no)
# var(mu2hat_no)
# mean(var2hat_no)