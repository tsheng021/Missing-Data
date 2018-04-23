library(MASS)
library(nlme)
library(reshape2)
library(parallel)
library(invgamma)
library(MCMCpack)
library(Hmisc)
library(gmodels)

rm(list = ls())


# mu_trt <- c(10, 20)
# mu_act <- c(10, 20)
# mu_pla <- c(10, 20)
mu_trt <- c(10, 20)
mu_act <- c(10, 10)
mu_pla <- c(10, 20)

vv <- matrix(c(1, .5, .5, 1),nrow=2)
N <- 100

# pi <- .7
# pi_c <- .15
# omega <- .7
# omega_c <- .15

pi <- .3
pi_c <- .3
omega <- .3
omega_c <- .3

# pi <- 1
# pi_c <- 0
# omega <- 1
# omega_c <- 0
# 
# pi <- .5
# pi_c <- 0
# omega <- .5
# omega_c <- 0

S <- 5000



no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl,varlist=c('S', 'pi_c', 'pi', 'omega', 'omega_c', 'N', 'vv', 'mu_trt', 'mu_act', 'mu_pla'))
clusterSetRNGStream(cl, 12345)
clusterEvalQ(cl, {
  library(MASS)
  library(nlme)
  library(reshape2)
  library(invgamma)
  library(MCMCpack)
  library(Hmisc)
  library(gmodels)
})


result <- parSapply(cl, 1:S, function(x){
  ## Estimand J2Z ##
  n1 <- sample(1:3,N,replace=TRUE,prob=c(pi,pi_c,1-pi-pi_c))  # partition size
  n1_trt <- sum(n1==1)
  n1_act <- sum(n1==2)
  n1_mis <- sum(n1==3)
  
  if(n1_trt*n1_act*n1_mis == 0){
    n1 <- sample(1:3,N,replace=TRUE,prob=c(pi,pi_c,1-pi-pi_c))  # partition size
    n1_trt <- sum(n1==1)
    n1_act <- sum(n1==2)
    n1_mis <- sum(n1==3)
  }
  
  
  y_trt <- y_act <- y_mis <- data.frame()
  
  ## treatment on treatment ##
  y_trt <- data.frame(id=as.character(1:n1_trt), 
                      mvrnorm(n1_trt, mu_trt, vv), 
                      dose = 'treatment', group = 'one')
  colnames(y_trt) <- c('id', 't1', 't2', 'dose', 'group')
  
  
  ## treatment on active control ##
  y_act <- data.frame(id=as.character((n1_trt+1):(n1_trt+n1_act)), 
                      mvrnorm(n1_act, mu_act, vv), 
                      dose = 'active', group='one')
  colnames(y_act) <- c('id', 't1', 't2', 'dose', 'group')
  
  
  
  
  ## treatment missing Case ##
  y_mis <- data.frame(id=as.character((n1_trt+n1_act+1):N),
                      mvrnorm(n1_mis, mu_act, vv), 
                      dose='active', group='one')
  colnames(y_mis) <- c('id', 't1', 't2', 'dose', 'group')  
  
  
  y_hyp <- rbind(y_trt, y_act, y_mis)  # The hypothetical complete data
  y_obs <- rbind(y_trt, y_act)   # The observed complete data
  
  
  r <- n1_trt
  q <- n1_trt + n1_act
  
  
  
  # If completely observed
  theta <- mean(y_hyp$t2)
  vartotal <- var(y_hyp$t2)
  se <- sqrt(vartotal/N)
  
  # If partial observed
  theta_obs <- mean(y_obs$t2)
  vartotal_obs <- var(y_obs$t2)
  se_obs <- sqrt(vartotal_obs/q)
  
  ## MMRM ##
  y_obs_long <- melt(y_obs)
  colnames(y_obs_long) <- c('id','dose','group','time','val')
  
  m1 <- lme(val~time*dose, 
            random=~1+time|id, 
            method='REML', 
            data=y_obs_long,
            control=lmeControl(returnObject=TRUE))
  
  
  
  ## MI Rubin ##
  MI <- 10
  z <- matrix(c(1,0,1,0,1,1,1,1), nrow=2, byrow=T)
  mu_act <- z%*%t(t(coef(summary(m1))[,1]))
  
  z <- matrix(c(1,1,0,1), nrow=2)
  sigma_mu <- z%*%getVarCov(m1)%*%t(z)/n1_act
  
  theta_mi_j <- vartotal_mi_j <- varcond_mi_j <- rep(NA, MI)
  for(j in 1:MI){
    mu1_tilde <- mean(y_mis$t1)
    mu_act_tilde <- c(mu1_tilde, 
                      rnorm(1, mu_act[2]+sigma_mu[1,2]/sigma_mu[1,1]*(mu1_tilde-mu_act[1]),
                            sqrt(sigma_mu[2,2]-sigma_mu[1,2]/sigma_mu[1,1]*sigma_mu[1,2])))
    
    
    sig_tilde2 <- rinvchisq(1, N-n1_mis-1)*m1$sigma^2*(N-n1_mis-1)
    
    y_mis$t2 <- rnorm(n1_mis,
                      mu_act_tilde[2],
                      sqrt(sig_tilde2))
    
    y2_mi <- c(y_obs$t2, y_mis$t2)
    theta_mi_j[j] <- mean(y2_mi)
    vartotal_mi_j[j] <- var(y2_mi)/N
  }
  theta_mi <- mean(theta_mi_j)
  W <- mean(vartotal_mi_j)
  B <- var(theta_mi_j)
  vartotal_mi <-  W + (1+1/MI)*B
  se_mi <- sqrt(vartotal_mi)
  
  r <- c(theta, vartotal, se,
         theta_obs, vartotal_obs, se_obs,
         theta_mi, vartotal_mi, se_mi)
})

stopCluster(cl)

rowMeans(result)
sqrt(apply(result,1,var) )
