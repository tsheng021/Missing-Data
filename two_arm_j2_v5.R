library(MASS)
library(nlme)
library(reshape2)
library(parallel)
library(invgamma)
library(MCMCpack)
library(Hmisc)
library(gmodels)

rm(list = ls())

mu_trt <- c(1, 4)
mu_act <- c(.5, 2)
mu_mis <- c(.2, 2)
eta_pla <- c(1, 20/3)
eta_act <- c(.5, 2)
eta_mis <- c(.2, 2)

vv <- matrix(c(4, 3, 3, 9),nrow=2)
N <- 100

# pi <- .7
# pi_c <- .15
# omega <- .7
# omega_c <- .15

# pi <- .3
# pi_c <- .3
# omega <- .3
# omega_c <- .3

pi <- .7
pi_c <- .15
omega <- .3
omega_c <- .3

S <- 5000

no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl,varlist=c('S', 'N', 'vv', 
                           'pi_c', 'pi', 'omega', 'omega_c', 
                           'mu_trt', 'mu_act', 'mu_mis',
                           'eta_pla', 'eta_act', 'eta_mis'))
clusterSetRNGStream(cl, 1234)
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
  
  while(n1_trt*n1_act*n1_mis == 0){
    n1 <- sample(1:3,N,replace=TRUE,prob=c(pi,pi_c,1-pi-pi_c))  # partition size
    n1_trt <- sum(n1==1)
    n1_act <- sum(n1==2)
    n1_mis <- sum(n1==3)
  }
  
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
                      mvrnorm(n1_mis, mu_mis, vv), 
                      dose='miss', group='one')
  colnames(y_mis) <- c('id', 't1', 't2', 'dose', 'group')  
  
  
  y_hyp <- rbind(y_trt, y_act, y_mis)  # The hypothetical complete data
  y_obs <- rbind(y_trt, y_act)   # The observed complete data
  
  
  r <- n1_trt
  q <- n1_trt + n1_act
  
  
  n2 <- sample(1:3,N,replace=TRUE,prob=c(omega,omega_c,1-omega-omega_c))  # partition size
  n2_pla <- sum(n2==1)
  n2_act <- sum(n2==2)
  n2_mis <- sum(n2==3)
  
  while(n2_pla*n2_act*n2_mis == 0){
    n2 <- sample(1:3,N,replace=TRUE,prob=c(omega,omega_c,1-omega-omega_c))  # partition size
    n2_pla <- sum(n2==1)
    n2_act <- sum(n2==2)
    n2_mis <- sum(n2==3)
  }
  
  z_pla <- z_act <- z_mis <- data.frame()
  ## placebo on placebo ##
  z_pla <- data.frame(id=as.character(1:n2_pla+N), 
                      mvrnorm(n2_pla, eta_pla, vv), 
                      dose = 'placebo', group = 'two')
  colnames(z_pla) <- c('id', 't1', 't2', 'dose', 'group')
  
  
  ## treatment on active control ##
  z_act <- data.frame(id=as.character((n2_pla+1):(n2_pla+n2_act)+N), 
                      mvrnorm(n2_act, eta_act, vv), 
                      dose = 'active', group='two')
  colnames(z_act) <- c('id', 't1', 't2', 'dose', 'group')
  
  
  
  
  ## treatment missing Case ##
  z_mis <- data.frame(id=as.character((n2_pla+n2_act+1):N+N),
                      mvrnorm(n2_mis, eta_mis, vv), 
                      dose='miss', group='two')
  colnames(z_mis) <- c('id', 't1', 't2', 'dose', 'group')  
  
  
  z_hyp <- rbind(z_pla, z_act, z_mis)  # The hypothetical complete data
  z_obs <- rbind(z_pla, z_act)   # The observed complete data
  
  
  p <- n2_pla
  s <- n2_pla + n2_act
  
  # If completely observed
  diff <- mean(y_hyp$t2)-mean(z_hyp$t2)
  vartotal <- var(y_hyp$t2)/N+var(z_hyp$t2)/N
  se <- sqrt(vartotal)
  # pval <- abs(diff/se)>1.96
  pval <- t.test(y_hyp$t2, z_hyp$t2, var.equal=T)$p.value<.05
  
  # If partial observed
  diff_obs <- mean(y_obs$t2)-mean(z_obs$t2)
  vartotal_obs_trt <- var(y_obs$t2)/q
  vartotal_obs_pla <- var(z_obs$t2)/s
  vartotal_obs <- vartotal_obs_trt + vartotal_obs_pla
  se_obs <- sqrt(vartotal_obs)
  # pval_obs <- abs(diff_obs/se_obs)>1.96
  pval_obs <- t.test(y_obs$t2, z_obs$t2, var.equal=T)$p.value<.05
  
  ## MMRM ##
  y_obs_long <- melt(y_obs)
  colnames(y_obs_long) <- c('id','dose','group','time','val')
  z_obs_long <- melt(z_obs)
  colnames(z_obs_long) <- c('id','dose','group','time','val')
  
  comp <- rbind(y_obs_long, z_obs_long)
  comp <- within(comp,  ind <- as.factor(paste(dose,group, sep='_')))
  
  m1 <- lme(val~time*ind, 
            random=~1+time|id, 
            method='REML', 
            data=comp,
            control=lmeControl(returnObject=TRUE))
  
  
  
  ## MI Rubin ##
  MI <- 50
  mu_act <- colMeans(y_act[,2:3])
  eta_act <- colMeans(z_act[,2:3])
  z <- matrix(c(1,1,0,1), nrow=2)
  sigma <- z%*%getVarCov(m1)%*%t(z) + diag(2)*m1$sigma^2
  
  diff_mi_j <- vartotal_mi_j <- rep(NA, MI)
  for(j in 1:MI){
    mu_act_tilde <- mvrnorm(1, mu_act, sigma/(n1_act-1))
    eta_act_tilde <- mvrnorm(1, eta_act, sigma/(n2_act-1))
    sigma_tilde <- riwish(s+q-6, sigma)*(s+q-6)
    
    y_mis$t2 <- rnorm(n1_mis,
                      mu_act_tilde[2]+sigma_tilde[1,2]/sigma_tilde[1,1]*(y_mis$t1-mu_act_tilde[1]),
                      sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
    z_mis$t2 <- rnorm(n2_mis,
                      eta_act_tilde[2]+sigma_tilde[1,2]/sigma_tilde[1,1]*(z_mis$t1-eta_act_tilde[1]),
                      sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
    
    y2_mi <- c(y_obs$t2, y_mis$t2)
    z2_mi <- c(z_obs$t2, z_mis$t2)
    diff_mi_j[j] <- mean(y2_mi)-mean(z2_mi)
    vartotal_mi_j[j] <- var(y2_mi)/N + var(z2_mi)/N
  }
  diff_mi <- mean(diff_mi_j)
  W <- mean(vartotal_mi_j)
  B <- var(diff_mi_j)
  vartotal_mi <-  W + (1+1/MI)*B
  se_mi <- sqrt(vartotal_mi)
  df <- (1+MI/(MI+1)*W/B)*(MI-1)
  t <- abs(diff_mi/se_mi)
  pval_mi <- (2*pt(-abs(t),df))<.05
  
  
  r <- c(diff, vartotal, se, pval, 
         diff_obs, vartotal_obs, se_obs, pval_obs,
         diff_mi, vartotal_mi, se_mi, pval_mi)
})

stopCluster(cl)

rowMeans(result)
sqrt(apply(result,1,var) )
