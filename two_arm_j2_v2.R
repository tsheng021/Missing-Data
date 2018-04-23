library(MASS)
library(condMVNorm)
library(nlme)
library(reshape2)
library(parallel)
library(invgamma)
library(MCMCpack)
library(Hmisc)
library(gmodels)

rm(list = ls())


mu_trt <- c(10, 20)
mu_act <- c(10, 20)
mu_pla <- c(10, 20)

vv <- matrix(c(1, .5, .5, 1),nrow=2)
N <- 400
pi <- .8
pi_c <- .15
omega <- .8
omega_c <- .15

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
  library(condMVNorm)
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
  n2 <- sample(1:3,N,replace=TRUE,prob=c(omega,omega_c,1-omega-omega_c))  # partition size
  n1_trt <- sum(n1==1)
  n1_act <- sum(n1==2)
  n1_mis <- sum(n1==3)
  n2_pla <- sum(n2==1)
  n2_act <- sum(n2==2)
  n2_mis <- sum(n2==3)
  
  y_trt <- y_act <- y_mis <- z_pla <- z_act <- z_mis <- data.frame()
  ## treatment on treatment ##
  if (n1_trt>0){
    y_trt <- data.frame(id=as.character(1:n1_trt), 
                        mvrnorm(n1_trt, mu_trt, vv), 
                        dose = 'treatment', group = 'one')
    colnames(y_trt) <- c('id', 't1', 't2', 'dose', 'group')
  }
  
  ## treatment on active control ##
  if (n1_act>0){
    y_act <- data.frame(id=as.character((n1_trt+1):(n1_trt+n1_act)), 
                        mvrnorm(n1_act, mu_act, vv), 
                        dose = 'active', group='one')
    colnames(y_act) <- c('id', 't1', 't2', 'dose', 'group')
  }
  
  
  ## treatment missing Case ##
  if(n1_mis>0){
    y_mis <- data.frame(id=as.character((n1_trt+n1_act+1):N),
                        mvrnorm(n1_mis, mu_act, vv), 
                        dose='active', group='one')
    colnames(y_mis) <- c('id', 't1', 't2', 'dose', 'group')
  }
  
  
  
  ## placebo on placebo ##
  if (n2_pla>0){
    z_pla <- data.frame(id=as.character(1:n2_pla+N), 
                        mvrnorm(n2_pla, mu_pla, vv), 
                        dose = 'placebo', group = 'two')
    colnames(z_pla) <- c('id', 't1', 't2', 'dose', 'group')
  }
  
  ## placebo on active control ##
  if (n2_act>0){
    z_act <- data.frame(id=as.character((n2_pla+1):(n2_pla+n2_act)+N), 
                        mvrnorm(n2_act, mu_act, vv), 
                        dose = 'active', group='two')
    colnames(z_act) <- c('id', 't1', 't2', 'dose', 'group')
  }
  
  ## placebo missing Case ##
  if(n2_mis>0){
    z_mis <- data.frame(id=as.character((n2_pla+n2_act+1):N+N),
                        mvrnorm(n2_mis, mu_act, vv), 
                        dose='active', group='two')
    colnames(z_mis) <- c('id', 't1', 't2', 'dose', 'group')
  }
  
  
  # The hypothetical complete data #
  y_hyp <- rbind(y_trt, y_act, y_mis)
  z_hyp <- rbind(z_pla, z_act, z_mis)
  # The observed complete data #
  y_obs <- rbind(y_trt, y_act)
  z_obs <- rbind(z_pla, z_act)
  
  r <- n1_trt
  q <- n1_trt + n1_act
  p <- n2_pla
  s <- n2_pla + n2_act
  
  
  # If completely observed
  diffhat <- mean(y_hyp$t2)-mean(z_hyp$t2)
  if (!is.null(c(y_act$t2, y_mis$t2, z_act$t2, z_mis$t2)))
    sig2hat <- ((r-1)*var(y_trt$t2) + (p-1)*var(z_pla$t2) + (N-r+N-p-1)*var(c(y_act$t2, y_mis$t2, z_act$t2, z_mis$t2)))/(2*N-3)
  if (is.null(c(y_act$t2, y_mis$t2, z_act$t2, z_mis$t2)))
    sig2hat <- ((r-1)*var(y_trt$t2) + (p-1)*var(z_pla$t2))/(2*N-2)
  vardiffhat <- var(y_hyp$t2)/N+var(z_hyp$t2)/N
  pval <- t.test(y_hyp$t2, z_hyp$t2)$p.value<.05
  
  
  ## MMRM ##
  y_obs_long <- melt(y_obs)
  colnames(y_obs_long) <- c('id','dose','group','time','val')
  z_obs_long <- melt(z_obs)
  colnames(z_obs_long) <- c('id','dose','group','time','val')
  comp <- rbind(y_obs_long, z_obs_long)
  
  z <- matrix(c(1,1,0,1), nrow=2)
  m <- lme(val~time*dose, 
           random=~1+time|id, 
           method='REML', 
           data=comp,
           control=lmeControl(returnObject=TRUE))
  
  if (q==r & p==s)
    cm <- rbind('contrast1'=c(0,0,-p/N,-p/N))
  if (q!=r | p!=s)
    cm <- rbind('contrast1'=c(0,0,(p-r)/N,-p/N,(p-r)/N,-p/N))
  cont <- estimable(m,cm)
  diffhat_mmrm <- cont[,1]
  vardiffhat_mmrm <- cont[,2]^2
  pval_mmrm <- cont[,5]<.05
  
  
  ## MI Rubin ##
  MI <- 20
  z <- matrix(c(1,1,0,1), nrow=2)
  Sigmahat <- z%*%getVarCov(m)%*%t(z)+diag(2)*m$sigma^2
  
  if (q==r & p==s){
    mu_trt_hat <- matrix(c(1,1,0,1,0,0,0,1), nrow=2)%*%t(t(coef(summary(m))[,1]))
    mu_pla_hat <- matrix(c(1,1,0,1,1,1,0,1), nrow=2)%*%t(t(coef(summary(m))[,1]))
  }
  if (q!=r | p!=s){
    mu_trt_hat <- matrix(c(1,0,0,0,0,0,1,1,0,0,0,0), nrow=2, byrow=T)%*%t(t(coef(summary(m))[,1]))
    mu_act_hat <- matrix(c(1,0,1,0,0,0,1,1,1,0,1,0), nrow=2, byrow=T)%*%t(t(coef(summary(m))[,1]))
    mu_pla_hat <- matrix(c(1,0,0,1,0,1,1,1,0,1,0,1), nrow=2, byrow=T)%*%t(t(coef(summary(m))[,1]))
  }
  n_trt <- r
  n_pla <- p
  n_act <- s-p+q-r
  
  diff_mi_j <- var_mi_j <- rep(NA, MI) 
  for(j in 1:MI){
    # 
    # mu_trt_tilde <- mvrnorm(1, mu_trt, Sigmahat/(n_trt-1))
    # mu_pla_tilde <- mvrnorm(1, mu_pla, Sigmahat/(n_pla-1))
    mu_act_tilde <- mvrnorm(1, mu_act_hat, Sigmahat/(n_act-1))
    
    
    sig_tilde <- rinvchisq(1, 2*N-3)*m$sigma^2*(2*N-3)
    
    Sigma_tilde <- riwish(2*N-3, Sigmahat)*(2*N-3)

    
    
    y_mis$t2 <- rnorm(n1_mis, 
                      mu_act_tilde[2]+Sigma_tilde[1,2]/Sigma_tilde[1,1]*(y_mis$t1-mu_act_tilde[1]), 
                      sqrt(Sigma_tilde[2,2]-Sigma_tilde[1,2]/Sigma_tilde[1,1]*Sigma_tilde[1,2]))
    z_mis$t2 <- rnorm(n2_mis, 
                      mu_act_tilde[2]+Sigma_tilde[1,2]/Sigma_tilde[1,1]*(z_mis$t1-mu_act_tilde[1]), 
                      sqrt(Sigma_tilde[2,2]-Sigma_tilde[1,2]/Sigma_tilde[1,1]*Sigma_tilde[1,2]))
    y2_mi <- c(y_obs$t2, y_mis$t2)
    z2_mi <- c(z_obs$t2, z_mis$t2)
    diff_mi_j[j] <- mean(y2_mi-z2_mi)
    var_mi_j[j] <- var(y2_mi)/N+var(z2_mi)/N
  }
  diffhat_mi <- mean(diff_mi_j)
  W <- mean(var_mi_j)
  B <- var(diff_mi_j)
  vardiffhat_mi <- W + (1+1/MI)*B
  df <- (1+MI/(MI+1)*W/B)*(MI-1)
  t <- abs(diffhat_mi/sqrt(vardiffhat_mi))
  pval_mi <- (2*pt(-abs(t),df))<.05
  
  r <- c(diffhat, sig2hat,vardiffhat, pval, 
         diffhat_mmrm, vardiffhat_mmrm, pval_mmrm,
         diffhat_mi, vardiffhat_mi,pval_mi)
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