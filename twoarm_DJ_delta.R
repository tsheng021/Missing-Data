start <- proc.time()
#setwd('./Documents/Missing Data/program/result/')
library(parallel)
library(MASS)
library(reshape2)
library(plyr)

# Simulation parameters
seed <- 12345
S <- 100
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

# MNAR4
# psi_p <- c(-4.0, 0, 0.4)
# psi_d <- c(-1.55, 0, -0.4)
# MAR1
# psi_p <- c(-2.6, -0.2, 0)
psi_cp <- c(-2.0, -0.2, 0)
psi_mp <- c(-2.0, -0.2, 0)
psi_cd <- c(-2.0, -0.2, 0)
psi_md <- c(-2.0, -0.2, 0)

# psi_p <- c(-4.0, 0, 0.4)
# psi_d <- c(-2.0, -0.2, 0)
true_dj=0


no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl,varlist=c('S', 'alpha','N','t',
                           'mu_pj', 'mu_dj', 'mu_pf', 'mu_df', 'vv',
                           'psi_cp', 'psi_mp','psi_cd','psi_md',
                           'true_dj',
                           'seed'))
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
  while (length(table(ind_p[,3]))<3){
    for(i in 1:nrow(y_p)){
      for(j in 3:t){
        pc_j <- 1/(1+exp(-(psi_cp[1]+psi_cp[2]*y_p[i,j-1]+psi_cp[3]*y_p[i,j])))
        pm_j <- 1/(1+exp(-(psi_mp[1]+psi_mp[2]*y_p[i,j-1]+psi_mp[3]*y_p[i,j])))
        ind_p[i,j-2] <- sample(1:3, 1, prob=c(1-pc_j, pc_j-pc_j*pm_j, pc_j*pm_j))
        if (j>3){
          if(ind_p[i,j-2]<ind_p[i,j-3])  ind_p[i,j-2] <- ind_p[i,j-3]
        } 
      }
    }
  }
  ind_p <- cbind(1,1,ind_p)
  
  ind_d <- matrix(rep(NA,(t-2)*N),nrow=N)
  while (length(table(ind_d[,3]))<3){
    for(i in 1:nrow(y_d)){
      for(j in 3:t){
        pc_j <- 1/(1+exp(-(psi_cd[1]+psi_cd[2]*y_d[i,j-1]+psi_cd[3]*y_d[i,j])))
        pm_j <- 1/(1+exp(-(psi_md[1]+psi_md[2]*y_d[i,j-1]+psi_md[3]*y_d[i,j])))
        ind_d[i,j-2] <- sample(1:3, 1, prob=c(1-pc_j, pc_j-pc_j*pm_j, pc_j*pm_j))
        if (j>3){
          if(ind_d[i,j-2]<ind_d[i,j-3])  ind_d[i,j-2] <- ind_d[i,j-3]
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
    obs <- which(ind_p[i,]==1)
    mis <- which(ind_p[i,]>1)
    if(length(mis)>0){
      mu <- mu_pf[mis] + vv[mis, obs]%*%solve(vv[obs, obs])%*%t(x[obs]-mu_pj[obs])
      v <- vv[mis, mis] - vv[mis, obs]%*%solve(vv[obs, obs])%*%vv[obs, mis]
      x[mis] <- mvrnorm(1, mu, v)
      pla_df[i,-1] <- x
    }
    switch_ind[i] <- max(which(ind_p[i,]==1))
    miss_ind[i] <- min(which(ind_p[i,]==3),5)
  }
  pla_df <- data.frame(pla_df, switch_ind, miss_ind)
  
  trt_df <- data.frame(1:N+N,y_d)
  colnames(trt_df) <- c('ID','bl','t1', 't2', 't3','t4')
  switch_ind <- miss_ind <- rep(NA, N) 
  for (i in 1:N){
    x<- trt_df[i,-1]
    obs <- which(ind_d[i,]==1)
    mis <- which(ind_d[i,]>1)
    if(length(mis)>0){
      mu <- mu_df[mis] + vv[mis, obs]%*%solve(vv[obs, obs])%*%t(x[obs]-mu_dj[obs])
      v <- vv[mis, mis] - vv[mis, obs]%*%solve(vv[obs, obs])%*%vv[obs, mis]
      x[mis] <- mvrnorm(1, mu, v)
      trt_df[i,-1] <- x
    }
    switch_ind[i] <- max(which(ind_d[i,]==1))
    miss_ind[i] <- min(which(ind_p[i,]==3),5)
  }
  trt_df <- data.frame(trt_df, switch_ind, miss_ind)
  
  # Change the data format from wide to long
  pla_df_long <- cbind.data.frame(melt(pla_df,id.var=c('ID', 'bl', 'switch_ind', 'miss_ind')), 
                                  dose='placebo',stringsAsFactors = FALSE)
  pla_df_long$dose[which(as.numeric(pla_df_long$variable)>=pla_df_long$switch_ind)] <- 'rescue'
  trt_df_long <- cbind.data.frame(melt(trt_df,id.var=c('ID', 'bl', 'switch_ind', 'miss_ind')), 
                                  dose='treatment',stringsAsFactors = FALSE)
  trt_df_long$dose[which(as.numeric(trt_df_long$variable)>=trt_df_long$switch_ind)] <- 'rescue'
  
  tot_df_long <- rbind(cbind(pla_df_long,arm=1),cbind(trt_df_long,arm=2))
  tot_df_long$delta <- as.numeric(tot_df_long$dose=='rescue')
  # `De Facto` estimator
  # Complete
  ## De Jure complete
  theta_dj <- mean(subset(pla_df_long,variable=='t4'&dose=='placebo')$value)-
    mean(subset(trt_df_long,variable=='t4'&dose=='treatment')$value)
  vartotal_dj <- var(subset(pla_df_long,variable=='t4'&dose=='placebo')$value)+
    var(subset(trt_df_long,variable=='t4'&dose=='treatment')$value)
  df <- (length(subset(pla_df_long,variable=='t4'&dose=='placebo')$value) + 
           length(subset(trt_df_long,variable=='t4'&dose=='treatment')$value))/2
  se_dj <- sqrt(vartotal_dj/df)
  low_dj <- theta_dj - qt(1-alpha/2, df-1)*se_dj
  upp_dj <- theta_dj + qt(1-alpha/2, df-1)*se_dj
  cover_dj <- low_dj<true_dj & true_dj<upp_dj
  
  out <- c(out, theta_dj, se_dj, cover_dj)
  # MMRM
  m1 <- lme(value~arm+variable+delta+arm:variable+arm:delta, 
            random=~1+variable|ID, 
            method='REML', 
            data=tot_df_long,
            control=lmeControl(returnObject=TRUE))  
  est <- estimable(m1,
                   cm=cbind(
                     'arm' = -1,
                     'arm:variablet4' = -1
                   ),
                   conf.int=0.95)
  est <- as.numeric(est)
  theta_dj_mmrm <- est[1]
  se_dj_mmrm <- est[2]
  low_dj_mmrm <- est[6]
  upp_dj_mmrm <- est[7]
  cover_dj_mmrm <- low_dj_mmrm<true_dj & true_dj<upp_dj_mmrm
  out <- c(out, theta_dj_mmrm, se_dj_mmrm, cover_dj_mmrm)
  
  ## MI
  
  MI <- 50
  z <- matrix(c(1,0,0,0,
                1,1,0,0,
                1,1,1,0,
                1,1,1,1), nrow=4, byrow=T)
  sigma <- z%*%getVarCov(m1)%*%t(z) + diag(4)*m1$sigma^2
  z_arm <- matrix(c(1,0,0,0,0,0,0,0,0,0,
                    1,0,1,0,0,0,1,0,0,0,
                    1,0,0,1,0,0,0,1,0,0,
                    1,0,0,0,1,0,0,0,1,0,
                    1,1,0,0,0,0,0,0,0,0,
                    1,1,1,0,0,0,1,0,0,0,
                    1,1,0,1,0,0,0,1,0,0,
                    1,1,0,0,1,0,0,0,1,0), nrow=8,byrow=T)
  sigma_arm <- z_arm%*%m1$varFix%*%t(z_arm)
  mu_arm <- z_arm %*% summary(m1)$coef$fixed
  theta_dj_mi_j <- vartotal_dj_mi_j <- rep(NA, MI)
  pla_mis <- pla_df
  trt_mis <- trt_df
  for(j in 1:MI){
    mu_tilde <- mvrnorm(1, mu_arm, sigma_arm)
    mu_pla_tilde <- mu_tilde[1:4]
    mu_trt_tilde <- mu_tilde[5:8]
    sigma_tilde <- riwish(N-1, sigma)*(N-1)
    for(i in 1:N){
      ind <- pla_mis$miss_ind[i]-1
      if(ind<4){
        obs <- 1:(ind-1)
        mis <- ind:4
        mu <- mu_pla_tilde[mis] + sigma_tilde[mis, obs]%*%solve(sigma_tilde[obs, obs])%*%t(pla_mis[i,obs+2]-mu_pla_tilde[obs])
        v <- vv[mis, mis] - vv[mis, obs]%*%solve(vv[obs, obs])%*%vv[obs, mis]
        pla_mis[i,mis+2] <- mvrnorm(1, mu, v)
      }
    }
    for(i in 1:N){
      ind <- trt_mis$miss_ind[i]-1
      if(ind<4){
        obs <- 1:(ind-1)
        mis <- ind:4
        mu <- mu_trt_tilde[mis] + sigma_tilde[mis, obs]%*%solve(sigma_tilde[obs, obs])%*%t(trt_mis[i,obs+2]-mu_pla_tilde[obs])
        v <- vv[mis, mis] - vv[mis, obs]%*%solve(vv[obs, obs])%*%vv[obs, mis]
        trt_mis[i,mis+2] <- mvrnorm(1, mu, v)
      }
    }
    
    theta_dj_mi_j[j] <- mean(trt_mis$t2) - mean(pla_mis$t2)
    vartotal_dj_mi_j[j] <- (var(trt_mis$t2) + var(pla_mis$t2))/(N-1)
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

end <- proc.time()
round(rowMeans(summary),3)
round(sqrt(apply(summary,1,var)),3)