start <- proc.time()
#setwd('./Documents/Missing Data/program/result/')
library(parallel)
library(MASS)
library(reshape2)
library(plyr)

# Simulation parameters
seed <- 1234
S <- 5000
N <- 100
t <- 5
alpha <- 0.05
mu_pj <- c(0, 1.0, 1.8, 2.5, 3)
mu_dj <- c(0, 1.0, 1.8, 2.5, 3)
mu_pf <- c(0, 1.0/2, 1.8/2, 2.5/2, 3/2)
mu_df <- c(0, 1.0/2, 1.8/2, 2.5/2, 3/2)
# mu_d <- c(0, 1.3, 2.3, 3.2, 4)
sd <- c(2.0, 1.8, 2.0, 2.1, 2.2)
corr <- matrix(c(1.0, 0.6, 0.3, 0.2, 0.1,
                 0.6, 1.0, 0.7, 0.5, 0.2,
                 0.3, 0.7, 1.0, 0.6, 0.4,
                 0.2, 0.5, 0.6, 1.0, 0.5,
                 0.1, 0.2, 0.4, 0.5, 1.0), nrow=5)
vv <- diag(sd) %*% corr %*% diag(sd)


# S1
# psi_cp <- c(-2, 0, 0)
# psi_mp <- c(1, 0, 0)
# psi_cd <- c(-2, 0, 0)
# psi_md <- c(1, 0, 0)
# S2
# psi_cp <- c(-1.8, -0.2, 0)
# psi_mp <- c(2.4, -0.2, 0)
# psi_cd <- c(-1.8, -0.2, 0)
# psi_md <- c(2.4, -0.2, 0)
# S3
# psi_cp <- c(-1.2, -0.4, 0)
# psi_mp <- c(1.5, -0.4, 0)
# psi_cd <- c(-1.2, -0.4, 0)
# psi_md <- c(1.5, -0.4, 0)
# S4
# psi_cp <- c(-1.5, 0, -0.2)
# psi_mp <- c(2.2, 0, -0.2)
# psi_cd <- c(-1.5, 0, -0.2)
# psi_md <- c(2.2, 0, -0.2)
# S5
# psi_cp <- c(-1, 0, -0.4)
# psi_mp <- c(2, 0, -0.4)
# psi_cd <- c(-1, 0, -0.4)
# psi_md <- c(2, 0, -0.4)
# S6
# psi_cp <- c(-1.8, -0.4, 0)
# psi_mp <- c(2.4, -0.4, 0)
# psi_cd <- c(-2, 0.4, 0)
# psi_md <- c(-1, 0.4, 0)
# S7
psi_cp <- c(-1.5, 0, -0.4)
psi_mp <- c(2.2, 0, -0.4)
psi_cd <- c(-3, 0, 0.4)
psi_md <- c(-1, 0, 0.4)

pla_mis_rate_k <- pla_evt_rate_k <- trt_mis_rate_k <- trt_evt_rate_k <- true_df_k <- rep(NA, 5000)
for(k in 1:10){
  
  y_p <- mvrnorm(S, mu_pj, vv)
  y_d <- mvrnorm(S, mu_dj, vv)
  
  ind_p <- matrix(rep(NA,(t-2)*S),nrow=S)
  for(i in 1:nrow(y_p)){
    for(j in 3:t){
      pc_j <- 1/(1+exp(-(psi_cp[1]+psi_cp[2]*y_p[i,j-1]+psi_cp[3]*y_p[i,j])))
      pm_j <- 1/(1+exp(-(psi_mp[1]+psi_mp[2]*y_p[i,j-1]+psi_mp[3]*y_p[i,j])))
      ind_p[i,j-2] <- sample(1:3, 1, prob=c(1-pc_j, pc_j-pc_j*pm_j, pc_j*pm_j))
      if (ind_p[i,j-2] == 2) {
        ind_p[i,(j-2):3] = 2
        break
      }
      if (ind_p[i,j-2] == 3) {
        ind_p[i,(j-2):3] = 3
        break
      }
    }
  }
  while (min(which(colMeans(ind_p==3)>0))>1 | min(which(colMeans(ind_p==2)>0))>1){
    for(i in 1:nrow(y_p)){
      for(j in 3:t){
        pc_j <- 1/(1+exp(-(psi_cp[1]+psi_cp[2]*y_p[i,j-1]+psi_cp[3]*y_p[i,j])))
        pm_j <- 1/(1+exp(-(psi_mp[1]+psi_mp[2]*y_p[i,j-1]+psi_mp[3]*y_p[i,j])))
        ind_p[i,j-2] <- sample(1:3, 1, prob=c(1-pc_j, pc_j-pc_j*pm_j, pc_j*pm_j))
        if (ind_p[i,j-2] == 2) {
          ind_p[i,(j-2):3] = 2
          break
        }
        if (ind_p[i,j-2] == 3) {
          ind_p[i,(j-2):3] = 3
          break
        }
      }
    }
  }
  ind_p <- cbind(1,1,ind_p)
  
  
  ind_d <- matrix(rep(NA,(t-2)*S),nrow=S)
  for(i in 1:nrow(y_d)){
    for(j in 3:t){
      pc_j <- 1/(1+exp(-(psi_cd[1]+psi_cd[2]*y_d[i,j-1]+psi_cd[3]*y_d[i,j])))
      pm_j <- 1/(1+exp(-(psi_md[1]+psi_md[2]*y_d[i,j-1]+psi_md[3]*y_d[i,j])))
      ind_d[i,j-2] <- sample(1:3, 1, prob=c(1-pc_j, pc_j-pc_j*pm_j, pc_j*pm_j))
      if (ind_d[i,j-2] == 2) {
        ind_d[i,(j-2):3] = 2
        break
      }
      if (ind_d[i,j-2] == 3) {
        ind_d[i,(j-2):3] = 3
        break
      }
    }
  }
  while (min(which(colMeans(ind_d==3)>0))>1 | min(which(colMeans(ind_d==2)>0))>1){
    for(i in 1:nrow(y_d)){
      for(j in 3:t){
        pc_j <- 1/(1+exp(-(psi_cd[1]+psi_cd[2]*y_d[i,j-1]+psi_cd[3]*y_d[i,j])))
        pm_j <- 1/(1+exp(-(psi_md[1]+psi_md[2]*y_d[i,j-1]+psi_md[3]*y_d[i,j])))
        ind_d[i,j-2] <- sample(1:3, 1, prob=c(1-pc_j, pc_j-pc_j*pm_j, pc_j*pm_j))
        if (ind_d[i,j-2] == 2) {
          ind_d[i,(j-2):3] = 2
          break
        }
        if (ind_d[i,j-2] == 3) {
          ind_d[i,(j-2):3] = 3
          break
        }
      }
    }
  }
  ind_d <- cbind(1,1,ind_d)
  
  # Simulate true value
  # Placebo
  pla_df <- data.frame(1:S,y_p)
  colnames(pla_df) <- c('ID','bl','t1', 't2', 't3','t4')
  switch_ind <- miss_ind <- rep(NA, S)
  for (i in 1:S){
    x<- pla_df[i,-1]
    obs <- which(ind_p[i,]==1)[-1]
    mis <- which(ind_p[i,]>1)
    if(length(mis)>0){
      mu <- mu_pf[mis] + vv[mis, obs]%*%solve(vv[obs, obs])%*%t(x[obs]-mu_pj[obs])
      v <- vv[mis, mis] - vv[mis, obs]%*%solve(vv[obs, obs])%*%vv[obs, mis]
      x[mis] <- mvrnorm(1, mu, v)
      pla_df[i,-1] <- x
    }
    switch_ind[i] <- max(which(ind_p[i,]==1))
    miss_ind[i] <- min(which(ind_p[i,]==3),6)-1
  }
  pla_df <- data.frame(pla_df, switch_ind, miss_ind)
  
  trt_df <- data.frame(1:S+S,y_d)
  colnames(trt_df) <- c('ID','bl','t1', 't2', 't3','t4')
  switch_ind <- miss_ind <- rep(NA, S)
  for (i in 1:S){
    x<- trt_df[i,-1]
    obs <- which(ind_d[i,]==1)[-1]
    mis <- which(ind_d[i,]>1)
    if(length(mis)>0){
      mu <- mu_df[mis] + vv[mis, obs]%*%solve(vv[obs, obs])%*%t(x[obs]-mu_dj[obs])
      v <- vv[mis, mis] - vv[mis, obs]%*%solve(vv[obs, obs])%*%vv[obs, mis]
      x[mis] <- mvrnorm(1, mu, v)
      trt_df[i,-1] <- x
    }
    switch_ind[i] <- max(which(ind_d[i,]==1))
    miss_ind[i] <- min(which(ind_d[i,]==3),6)-1
  }
  trt_df <- data.frame(trt_df, switch_ind, miss_ind)
  
  pla_mis_rate_k <- mean(pla_df$miss_ind<5)
  pla_evt_rate_k <- mean(pla_df$switch_ind<5)
  
  trt_mis_rate_k <- mean(trt_df$miss_ind<5)
  trt_evt_rate_k <- mean(trt_df$switch_ind<5)
  
  
  true_df_k <- mean(trt_df$t4)-mean(pla_df$t4)
}

(pla_mis_rate <- mean(pla_mis_rate_k))
(pla_evt_rate <- mean(pla_evt_rate_k))
(trt_mis_rate <- mean(trt_mis_rate_k))
(trt_evt_rate <- mean(trt_evt_rate_k))
(true_df <- mean(true_df_k))
true_df
true_dj <- 0



S <- 1000
no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl,varlist=c('S', 'alpha','N','t', 'seed',
                           'mu_pj', 'mu_dj', 'mu_pf', 'mu_df', 'vv',
                           'psi_cp', 'psi_mp','psi_cd','psi_md',
                           'true_dj','true_df'
))
clusterEvalQ(cl, {
  library(MASS)
  library(reshape2)
  library(nlme)
  library(parallel)
  library(invgamma)
  library(MCMCpack)
  library(Hmisc)
  library(gmodels)
})
clusterSetRNGStream(cl, seed)

summary <- parSapply(cl, 1:S,  function(x){
  out <- c()
  # Simulate `De Jure' complete 
  y_p <- mvrnorm(N, mu_pj, vv)
  y_d <- mvrnorm(N, mu_dj, vv)
  
  ind_p <- matrix(rep(NA,(t-2)*N),nrow=N)
  for(i in 1:nrow(y_p)){
    for(j in 3:t){
      pc_j <- 1/(1+exp(-(psi_cp[1]+psi_cp[2]*y_p[i,j-1]+psi_cp[3]*y_p[i,j])))
      pm_j <- 1/(1+exp(-(psi_mp[1]+psi_mp[2]*y_p[i,j-1]+psi_mp[3]*y_p[i,j])))
      ind_p[i,j-2] <- sample(1:3, 1, prob=c(1-pc_j, pc_j-pc_j*pm_j, pc_j*pm_j))
      if (ind_p[i,j-2] == 2) {
        ind_p[i,(j-2):3] = 2
        break
      }
      if (ind_p[i,j-2] == 3) {
        ind_p[i,(j-2):3] = 3
        break
      }
    }
  }
  while (min(which(colMeans(ind_p==3)>0))>1 | min(which(colMeans(ind_p==2)>0))>1){
    for(i in 1:nrow(y_p)){
      for(j in 3:t){
        pc_j <- 1/(1+exp(-(psi_cp[1]+psi_cp[2]*y_p[i,j-1]+psi_cp[3]*y_p[i,j])))
        pm_j <- 1/(1+exp(-(psi_mp[1]+psi_mp[2]*y_p[i,j-1]+psi_mp[3]*y_p[i,j])))
        ind_p[i,j-2] <- sample(1:3, 1, prob=c(1-pc_j, pc_j-pc_j*pm_j, pc_j*pm_j))
        if (ind_p[i,j-2] == 2) {
          ind_p[i,(j-2):3] = 2
          break
        }
        if (ind_p[i,j-2] == 3) {
          ind_p[i,(j-2):3] = 3
          break
        }
      }
    }
  }
  ind_p <- cbind(1,1,ind_p)
  
  ind_d <- matrix(rep(NA,(t-2)*N),nrow=N)
  for(i in 1:nrow(y_d)){
    for(j in 3:t){
      pc_j <- 1/(1+exp(-(psi_cd[1]+psi_cd[2]*y_d[i,j-1]+psi_cd[3]*y_d[i,j])))
      pm_j <- 1/(1+exp(-(psi_md[1]+psi_md[2]*y_d[i,j-1]+psi_md[3]*y_d[i,j])))
      ind_d[i,j-2] <- sample(1:3, 1, prob=c(1-pc_j, pc_j-pc_j*pm_j, pc_j*pm_j))
      if (ind_d[i,j-2] == 2) {
        ind_d[i,(j-2):3] = 2
        break
      }
      if (ind_d[i,j-2] == 3) {
        ind_d[i,(j-2):3] = 3
        break
      }
    }
  }
  while (min(which(colMeans(ind_d==3)>0))>1 | min(which(colMeans(ind_d==2)>0))>1){
    for(i in 1:nrow(y_d)){
      for(j in 3:t){
        pc_j <- 1/(1+exp(-(psi_cd[1]+psi_cd[2]*y_d[i,j-1]+psi_cd[3]*y_d[i,j])))
        pm_j <- 1/(1+exp(-(psi_md[1]+psi_md[2]*y_d[i,j-1]+psi_md[3]*y_d[i,j])))
        ind_d[i,j-2] <- sample(1:3, 1, prob=c(1-pc_j, pc_j-pc_j*pm_j, pc_j*pm_j))
        if (ind_d[i,j-2] == 2) {
          ind_d[i,(j-2):3] = 2
          break
        }
        if (ind_d[i,j-2] == 3) {
          ind_d[i,(j-2):3] = 3
          break
        }
      }
    }
  }
  ind_d <- cbind(1,1,ind_d)
  
  # Simulate `De Facto` complete
  # Placebo 
  pla_df <- data.frame(1:N,y_p)
  colnames(pla_df) <- c('ID','bl','t1', 't2', 't3','t4')
  switch_ind <- miss_ind <- rep(NA, N) 
  for (i in 1:N){
    x<- pla_df[i,-1]
    obs <- which(ind_p[i,]==1)[-1]
    mis <- which(ind_p[i,]>1)
    if(length(mis)>0){
      mu <- mu_pf[mis] + vv[mis, obs]%*%solve(vv[obs, obs])%*%t(x[obs]-mu_pj[obs])
      v <- vv[mis, mis] - vv[mis, obs]%*%solve(vv[obs, obs])%*%vv[obs, mis]
      x[mis] <- mvrnorm(1, mu, v)
      pla_df[i,-1] <- x
    }
    switch_ind[i] <- max(which(ind_p[i,]==1))
    miss_ind[i] <- min(which(ind_p[i,]==3),6)-1
  }
  pla_df <- data.frame(pla_df, switch_ind, miss_ind)
  
  trt_df <- data.frame(1:N+N,y_d)
  colnames(trt_df) <- c('ID','bl','t1', 't2', 't3','t4')
  switch_ind <- miss_ind <- rep(NA, N) 
  for (i in 1:N){
    x<- trt_df[i,-1]
    obs <- which(ind_d[i,]==1)[-1]
    mis <- which(ind_d[i,]>1)
    if(length(mis)>0){
      mu <- mu_df[mis] + vv[mis, obs]%*%solve(vv[obs, obs])%*%t(x[obs]-mu_dj[obs])
      v <- vv[mis, mis] - vv[mis, obs]%*%solve(vv[obs, obs])%*%vv[obs, mis]
      x[mis] <- mvrnorm(1, mu, v)
      trt_df[i,-1] <- x
    }
    switch_ind[i] <- max(which(ind_d[i,]==1))
    miss_ind[i] <- min(which(ind_d[i,]==3),6)-1
  }
  trt_df <- data.frame(trt_df, switch_ind, miss_ind)
  
  # Change the data format from wide to long
  pla_df_long <- cbind.data.frame(melt(pla_df,id.var=c('ID', 'bl', 'switch_ind', 'miss_ind')), 
                                  dose='placebo',stringsAsFactors = FALSE)
  pla_df_long$dose[which(as.numeric(pla_df_long$variable)>=pla_df_long$switch_ind)] <- 'rescue'
  trt_df_long <- cbind.data.frame(melt(trt_df,id.var=c('ID', 'bl', 'switch_ind', 'miss_ind')), 
                                  dose='treatment',stringsAsFactors = FALSE)
  trt_df_long$dose[which(as.numeric(trt_df_long$variable)>=trt_df_long$switch_ind)] <- 'rescue'
  
  tot_df_long <- rbind(cbind(pla_df_long,arm=0),cbind(trt_df_long,arm=1))
  tot_df_long$delta <- as.numeric(tot_df_long$dose=='rescue')
  
  
  # `De Facto` estimator
  # Complete
  ## De Facto observed
  pla_obs <- pla_df$t4[pla_df$miss_ind==5]
  trt_obs <- trt_df$t4[trt_df$miss_ind==5]
  n1 <- length(pla_obs)
  n2 <- length(trt_obs)
  theta_df <- mean(trt_obs)-mean(pla_obs)
  vartotal_df <- ((n2-1)*var(trt_obs)+(n1-1)*var(pla_obs))/(n1+n2-2)
  se_df <- sqrt(vartotal_df*(1/n1+1/n2))
  low_df <- theta_df - qt(1-alpha/2, n1+n2-2)*se_df
  upp_df <- theta_df + qt(1-alpha/2, n1+n2-2)*se_df
  cover_df <- low_df<true_df & true_df<upp_df
  t <- abs(theta_df)/se_df
  df <- n1+n2-2
  pval_df <- ((2*pt(-t, df))<.05)*100
  out <- c(out, theta_df, se_df, cover_df, pval_df)
  # MMRM
  tot_mis_long <- tot_df_long
  tot_mis_long <- tot_mis_long[-which(tot_mis_long$miss_ind <= as.numeric(tot_mis_long$variable)),]
  tot_mis_long$d2t <- as.numeric(tot_mis_long$dose=='rescue' & tot_mis_long$variable=='t2' & tot_mis_long$arm == 1)
  tot_mis_long$d3t <- as.numeric(tot_mis_long$dose=='rescue' & tot_mis_long$variable=='t3' & tot_mis_long$arm == 1)
  tot_mis_long$d4t <- as.numeric(tot_mis_long$dose=='rescue' & tot_mis_long$variable=='t4' & tot_mis_long$arm == 1)
  tot_mis_long$d2p <- as.numeric(tot_mis_long$dose=='rescue' & tot_mis_long$variable=='t2' & tot_mis_long$arm == 0)
  tot_mis_long$d3p <- as.numeric(tot_mis_long$dose=='rescue' & tot_mis_long$variable=='t3' & tot_mis_long$arm == 0)
  tot_mis_long$d4p <- as.numeric(tot_mis_long$dose=='rescue' & tot_mis_long$variable=='t4' & tot_mis_long$arm == 0)
  cand <- c()
  if(sum(tot_mis_long$d2p)>0)
    cand <- c(cand, 'd2p')
  if(sum(tot_mis_long$d3p)>0)
    cand <- c(cand, 'd3p')
  if(sum(tot_mis_long$d4p)>0)
    cand <- c(cand, 'd4p')
  if(sum(tot_mis_long$d2t)>0)
    cand <- c(cand, 'd2t')
  if(sum(tot_mis_long$d3t)>0)
    cand <- c(cand, 'd3t')
  if(sum(tot_mis_long$d4t)>0)
    cand <- c(cand, 'd4t')
  m1 <- lme(as.formula(paste('value~arm+variable+arm:variable+', paste(cand, collapse='+'))),
            random=~1+variable|ID,
            method='REML',
            data=tot_mis_long,
            control=lmeControl(returnObject=TRUE))
  est <- estimable(m1,
                   cm=cbind(
                     'arm' = 1,
                     'arm:variablet4' = 1,
                     'd4t' =  mean(trt_df$switch_ind<5),
                     'd4p' = - mean(pla_df$switch_ind<5)
                   ),
                   conf.int=0.95)
  est <- as.numeric(est)
  theta_df_mmrm <- est[1]
  se_df_mmrm <- est[2]
  low_df_mmrm <- est[6]
  upp_df_mmrm <- est[7]
  cover_df_mmrm <- low_df_mmrm<true_df & true_df<upp_df_mmrm
  pval_df_mmrm <- (est[5]<.05)*100
  out <- c(out, theta_df_mmrm, se_df_mmrm, cover_df_mmrm, pval_df_mmrm)
  
  MI <- 50
  z <- matrix(c(1,0,0,0,
                1,1,0,0,
                1,0,1,0,
                1,0,0,1), nrow=4, byrow=T)
  sigma <- z%*%getVarCov(m1)%*%t(z) + diag(4)*m1$sigma^2

  cand_seg <- rep(0, length(cand))

  
  z_arm_djdf <-  matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                          1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                          1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                          1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                          1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                          1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                          1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0,
                          1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0,
                          1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1), nrow=14, byrow=T)

  sigma_arm_djdf <- z_arm_djdf%*%m1$varFix%*%t(z_arm_djdf)
  mu_arm_djdf <- z_arm_djdf %*% summary(m1)$coef$fixed
  theta_df_mi_j <- vartotal_df_mi_j <- rep(NA, MI)
  pla_mis <- pla_df
  trt_mis <- trt_df
  for(j in 1:MI){
    mu_djdf_tilde <- mvrnorm(1, mu_arm_djdf, sigma_arm_djdf)
    mu_pla_dj_tilde <-c(mu_djdf_tilde[1], mu_djdf_tilde[3:5])
    mu_trt_dj_tilde <- c(mu_djdf_tilde[2], mu_djdf_tilde[6:8])
    mu_pla_df_tilde <-c(mu_djdf_tilde[1], mu_djdf_tilde[9:11])
    mu_trt_df_tilde <- c(mu_djdf_tilde[2], mu_djdf_tilde[12:14])
    sigma_tilde <- riwish(2*N-2, sigma)*(2*N-2)
    for(i in 1:N){
      ind <- pla_mis$miss_ind[i]
      if(ind<5){
        obs <- (2:ind)-1
        mis <- ((ind+1):5)-1
        mu <- mu_pla_df_tilde[mis] + sigma_tilde[mis, obs]%*%solve(sigma_tilde[obs, obs])%*%t(pla_mis[i,obs+2]-mu_pla_dj_tilde[obs])
        v <- sigma_tilde[mis, mis] - sigma_tilde[mis, obs]%*%solve(sigma_tilde[obs, obs])%*%sigma_tilde[obs, mis]
        pla_mis[i,mis+2] <- mvrnorm(1, mu, v)
      }
    }
    for(i in 1:N){
      ind <- trt_mis$miss_ind[i]
      if(ind<5){
        obs <- (2:ind)-1
        mis <- ((ind+1):5)-1
        mu <- mu_trt_df_tilde[mis] + sigma_tilde[mis, obs]%*%solve(sigma_tilde[obs, obs])%*%t(trt_mis[i,obs+2]-mu_trt_dj_tilde[obs])
        v <- sigma_tilde[mis, mis] - sigma_tilde[mis, obs]%*%solve(sigma_tilde[obs, obs])%*%sigma_tilde[obs, mis]
        trt_mis[i,mis+2] <- mvrnorm(1, mu, v)
      }
    }

    theta_df_mi_j[j] <- mean(trt_mis$t4) - mean(pla_mis$t4)
    vartotal_df_mi_j[j] <- (var(trt_mis$t4) + var(pla_mis$t4))/(N-1)
  }
  theta_df_mi <- mean(theta_df_mi_j)
  W <- mean(vartotal_df_mi_j)
  B <- var(theta_df_mi_j)
  se_df_mi <-  sqrt(W + (1+1/MI)*B)
  df <- (1+MI/(MI+1)*W/B)*(MI-1)
  low_df_mi <- theta_df_mi - qt(1-alpha/2, df)*se_df_mi
  upp_df_mi <- theta_df_mi + qt(1-alpha/2, df)*se_df_mi
  cover_df_mi <- low_df_mi<true_df & true_df<upp_df_mi
  t <- abs(theta_df_mi/se_df_mi)
  pval_df_mi <- (2*pt(-t, df)<.05)*100
  out <- c(out, theta_df_mi, se_df_mi, cover_df_mi, pval_df_mi)
})

stopCluster(cl)

end <- proc.time()
round(rowMeans(summary,na.rm=T),2)
round(sqrt(apply(summary,1,function(x){var(x,na.rm=T)})),2)
