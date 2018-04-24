library(MASS)
library(nlme)
library(reshape2)
library(parallel)

rm(list = ls())


mu_dj <- c(0, 1, 4)
mu_df <- c(0, 1, 2)
mu_pj <- c(0, 1, 4)
mu_pf <- c(0, 1, 2)

sd <- c(2.0, 1.8, 2.2)
corr <- matrix(c(1.0, 0.6, 0.1,
                 0.6, 1.0, 0.2,
                 0.1, 0.2, 1.0), nrow=3)
vv <- diag(sd) %*% corr %*% diag(sd)
N <- 100
alpha <- .05
# S1
# psi1 <- c(.4, -.4, 0)
# psi2 <- c(-.3, -.4, 0)
# psi1p <- c(-.4, .4, 0)
# psi2p <- c(-1.6, .4, 0)
# S2
# psi1 <- c(.4, -.4, 0)
# psi2 <- c(.3, 0, -.4)
# psi1p <- c(-.4, .4, 0)
# psi2p  <- c(-1.8, 0, .4)
# S3
psi1 <- c(1.6, 0, -.4)
psi2 <- c(-.1, -.4, 0)
psi1p <- c(1.6, 0, -.4)
psi2p  <- c(-0.6, -.4, 0)
# S4
# psi1 <- c(1.6, 0, -.4)
# psi2 <- c(.3, 0, -.4)
# psi1p <- c(1.6, 0, -.4)
# psi2p  <- c(-1.8, 0, .4)

S <- 1000

df <- 
  switch_trt <- miss_trt <- 
  switch_pla <- miss_pla <- 
  c(NA,S)
for(i in 1:S){
  # Treatment arm
  y_dj <- mvrnorm(N, mu_dj, vv)
  y_dj <- data.frame(id=1:N, t0=y_dj[,1], t1=y_dj[,2], t2=y_dj[,3], arm='trt',stringsAsFactors=FALSE)
  px <- as.matrix(cbind(1,y_dj[,3:4])) %*% t(t(psi1))
  pt <- 1/(1+exp(-px))
  ind <- rbinom(N, 1, pt)+1
  while (sum(ind==2)==0) ind <- rbinom(N, 1, pt)+1
  ind_df <- which(ind==2)
  y_df <- y_dj
  y_df$arm <- 'trt'
  y_df$arm[ind_df] <- 'act'
  y_df$t2[ind_df] <- rnorm(length(ind_df), 
                           mean=mu_df[3]+vv[3,1:2]%*%solve(vv[1:2,1:2])%*%t(sweep(y_df[ind_df,2:3], 2, mu_dj[1:2])),
                           sd=sqrt(vv[3,3]-vv[3,1:2]%*%solve(vv[1:2,1:2])%*%vv[1:2,3]))
  switch_trt[i] <- length(ind_df)/N
  px <- as.matrix(cbind(1,y_df[ind_df,3:4])) %*% t(t(psi2))
  pt <- 1/(1+exp(-px))
  ind[ind==2] <- ind[ind==2] + rbinom(sum(ind==2), 1, pt)
  while (sum(ind==3)==0|sum(ind==2)==0) {
    ind[ind==3] <- ind[ind==3] - 1
    ind[ind==2]<-ind[ind==2] + rbinom(sum(ind==2), 1, pt)
  }
  
  ind_mis <- which(ind==3)
  miss_trt[i] <- length(ind_mis)/length(ind_df)
  
  
  # Placebo arm
  y_pj <- mvrnorm(N, mu_pj, vv)
  y_pj <- data.frame(id=1:N+N, t0=y_pj[,1], t1=y_pj[,2], t2=y_pj[,3], arm='pla',stringsAsFactors=FALSE)
  px <- as.matrix(cbind(1,y_pj[,3:4])) %*% t(t(psi1p))
  pt <- 1/(1+exp(-px))
  indp <- rbinom(N, 1, pt)+1
  while (sum(indp==2)==0) indp <- rbinom(N, 1, pt)+1
  y_pf <- y_pj
  indp_df <- which(indp==2)
  y_pf$arm <- 'pla'
  y_pf$arm[indp_df] <- 'act'
  y_pf$t2[indp_df] <- rnorm(length(indp_df), 
                            mean=mu_pf[3]+vv[3,1:2]%*%solve(vv[1:2,1:2])%*%t(sweep(y_pf[indp_df,2:3], 2, mu_pj[1:2])),
                            sd=sqrt(vv[3,3]-vv[3,1:2]%*%solve(vv[1:2,1:2])%*%vv[1:2,3]))
  switch_pla[i] <- length(indp_df)/N
  px <- as.matrix(cbind(1,y_pf[indp_df,3:4])) %*% t(t(psi2p))
  pt <- 1/(1+exp(-px))
  indp[indp==2] <- indp[indp==2] + rbinom(sum(indp==2), 1, pt)
  while (sum(indp==3)==0|sum(indp==2)==0) {
    indp[indp==3] <- ind[indp==3]-1
    indp[indp==2]<-indp[indp==2] + rbinom(sum(indp==2), 1, pt)
  }
  indp_mis <- which(indp==3)
  miss_pla[i] <- length(indp_mis)/length(indp_df)
  
  df[i] <- mean(y_df$t2)-mean(y_pf$t2)
}
true_df <- mean(df)
true_switch_trt <- mean(switch_trt)
true_miss_trt <- mean(miss_trt)
true_switch_pla <- mean(switch_pla)
true_miss_pla <- mean(miss_pla)
true_df
true_switch_trt
true_miss_trt
true_switch_pla
true_miss_pla


no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl,varlist=c('S', 'N', 'vv', 'alpha', 
                           'true_df',
                           'psi1', 'psi2', 
                           'psi1p', 'psi2p',
                           'mu_pj','mu_pf',
                           'mu_dj','mu_df'))
clusterSetRNGStream(cl, 31415)
clusterEvalQ(cl, {
  library(MASS)
  library(nlme)
  library(reshape2)
  library(mice)
  library(emmeans)
})


result <- parSapply(cl, 1:S, function(x){
  out <- c()
  y_dj <- mvrnorm(N, mu_dj, vv)
  y_dj <- data.frame(id=1:N, t0=y_dj[,1], t1=y_dj[,2], t2=y_dj[,3], arm='trt',stringsAsFactors=FALSE)
  px <- as.matrix(cbind(1,y_dj[,3:4])) %*% t(t(psi1))
  pt <- 1/(1+exp(-px))
  ind <- rbinom(N, 1, pt)+1
  while (sum(ind==2)==0) ind <- rbinom(N, 1, pt)+1
  ind_df <- which(ind==2)
  y_df <- y_dj
  y_df$arm[ind_df] <- 'tact'
  y_df$t2[ind_df] <- rnorm(length(ind_df), 
                           mean=mu_df[3]+vv[3,1:2]%*%solve(vv[1:2,1:2])%*%t(sweep(y_df[ind_df,2:3], 2, mu_dj[1:2])),
                           sd=sqrt(vv[3,3]-vv[3,1:2]%*%solve(vv[1:2,1:2])%*%vv[1:2,3]))
  px <- as.matrix(cbind(1,y_df[ind_df,3:4])) %*% t(t(psi2))
  pt <- 1/(1+exp(-px))
  ind[ind==2] <- 2 + rbinom(sum(ind==2), 1, pt)
  while (sum(ind==3)==0|sum(ind==2)==0){
    ind[ind==3] <- 2
    ind[ind==2]<- 2 + rbinom(sum(ind==2), 1, pt)
  } 
  ind_res <- which(ind==2)
  ind_mis <- which(ind==3)
  y_mis <- y_df
  y_mis$t2[ind_mis] <- NA
  y_mis_long <- melt(y_mis,id.vars = c('id','arm','t0'))
  colnames(y_mis_long) <- c('id','arm','bl','time','val')
  y_mis_long <- na.omit(y_mis_long)
  
  r <- N-length(ind_df)
  q <- N-r-length(ind_mis)
  
  
  # Placebo arm
  y_pj <- mvrnorm(N, mu_pj, vv)
  y_pj <- data.frame(id=1:N+N, t0=y_pj[,1], t1=y_pj[,2], t2=y_pj[,3], arm='pla',stringsAsFactors=FALSE)
  px <- as.matrix(cbind(1,y_pj[,3:4])) %*% t(t(psi1p))
  pt <- 1/(1+exp(-px))
  indp <- rbinom(N, 1, pt)+1
  while (sum(indp==2)==0) indp <- rbinom(N, 1, pt)+1
  y_pf <- y_pj
  indp_df <- which(indp==2)
  y_pf$arm[indp_df] <- 'pact'
  y_pf$t2[indp_df] <- rnorm(length(indp_df), 
                            mean=mu_pf[3]+vv[3,1:2]%*%solve(vv[1:2,1:2])%*%t(sweep(y_pf[indp_df,2:3], 2, mu_pj[1:2])),
                            sd=sqrt(vv[3,3]-vv[3,1:2]%*%solve(vv[1:2,1:2])%*%vv[1:2,3]))
  px <- as.matrix(cbind(1,y_pf[indp_df,3:4])) %*% t(t(psi2p))
  pt <- 1/(1+exp(-px))
  indp[indp==2] <- 2 + rbinom(sum(indp==2), 1, pt)
  while (sum(indp==3)==0|sum(indp==2)==0) {
    indp[indp==3] <- 2
    indp[indp==2] <- 2 + rbinom(sum(indp==2), 1, pt)
  } 
  indp_res <- which(indp==2)
  indp_mis <- which(indp==3)
  y_misp <- y_pf
  y_misp$t2[indp_mis] <- NA
  y_misp_long <- melt(y_misp,id.vars = c('id','arm','t0'))
  colnames(y_misp_long) <- c('id','arm','bl','time','val')
  y_misp_long <- na.omit(y_misp_long)
  
  l <- N-length(indp_df)
  m <- N-l-length(indp_mis)
  
  # Total dataset
  tot_mis_long <- rbind(data.frame(y_mis_long,group='treatment'), data.frame(y_misp_long,group='placebo'))
  
  # MMRM
  m1 <- lme(val~bl+time+arm+arm:time+bl:time,
            random=~(-1)+time|id,
            method='REML',
            data=tot_mis_long,
            control=lmeControl(returnObject=TRUE))
  
  
  bl_T <- mean(y_mis$t0)
  bl_P <- mean(y_misp$t0)
  
  # MAR-MMRM
  coeff <- summary(m1)$coef$fixed
  D_mar <- c(0, bl_T-bl_P, 0, -l/N, 1-r/N, r/N, -l/N, 1-r/N, r/N, bl_T-bl_P)
  theta_mar_mmrm <- D_mar%*%coeff
  se_mmrm <- sqrt(D_mar%*%vcov(m1)%*%D_mar)
  se_mar_mmrm <- sqrt(se_mmrm^2 +
                        r*(N-r)/N^3*((coeff[6]+coeff[9])^2+(coeff[5]+coeff[8])^2-2*(coeff[6]+coeff[9])*(coeff[5]+coeff[8]))+
                        l*(N-l)/N^3*(coeff[4]+coeff[7])^2)
  low_mar_mmrm <- theta_mar_mmrm - qt(1-alpha/2, N-2)*se_mar_mmrm
  upp_mar_mmrm <- theta_mar_mmrm + qt(1-alpha/2, N-2)*se_mar_mmrm
  cover_mar_mmrm <- low_mar_mmrm<true_df & true_df<upp_mar_mmrm
  t <- theta_mar_mmrm/se_mar_mmrm
  pval_mar_mmrm <- ((1-pt(t,2*N))<.025)*100
  out <- c(out, theta_mar_mmrm, se_mar_mmrm, cover_mar_mmrm, pval_mar_mmrm)
  
  # MI
  MI <- 20
  sigma <- getVarCov(m1) + diag(2)*m1$sigma^2
  theta_mar_mi_j <- vartotal_mar_mi_j <- 
    theta_mild_mi_j <- vartotal_mild_mi_j <- 
    theta_moderate_mi_j <- vartotal_moderate_mi_j <- 
    theta_extreme_mi_j <- vartotal_extreme_mi_j <-rep(NA, MI)
  for(j in 1:MI){
    sigma_tilde <- matrix(rWishart(1, 2*N-4, sigma),nrow=2)/(2*N-4)
    # MAR
    y_mis$t2[ind_mis] <- NA
    y_misp$t2[indp_mis] <- NA
    for(k in 1:N){
      if(is.na(y_mis$t2[k])){
        z_act <- matrix(c(1,y_mis$t0[k],0,0,1,0,0,0,0,0,
                          1,y_mis$t0[k],1,0,1,0,0,1,0,y_mis$t0[k]),nrow=2,byrow=T)
        mu_act <- z_act %*% coeff
        mu_act_tilde <- mvrnorm(1, mu_act, sigma)
        y_mis$t2[k] <- rnorm(1,
                             mu_act_tilde[2]+sigma_tilde[2,1]/sigma_tilde[1,1]*(y_mis$t1[k]-mu_act_tilde[1]),
                             sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
      }
      if(is.na(y_misp$t2[k])){
        z_act <- matrix(c(1,y_misp$t0[k],0,0,0,0,0,0,0,0,
                          1,y_misp$t0[k],1,0,0,0,0,0,0,y_misp$t0[k]),nrow=2,byrow=T)
        mu_act <- z_act %*% coeff
        mu_act_tilde <- mvrnorm(1, mu_act, sigma)
        y_misp$t2[k] <- rnorm(1,
                              mu_act_tilde[2]+sigma_tilde[2,1]/sigma_tilde[1,1]*(y_misp$t1[k]-mu_act_tilde[1]),
                              sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
        
      }
    }
    y_mi_tot <- rbind(y_mis,y_misp)
    y_mi <- melt(y_mi_tot,id.vars = c('id','t0','arm'))
    colnames(y_mi) <- c('id','bl','arm','time','val')
    lmi <- lme(val~bl+time+arm+arm:time+bl:time,
               random=~(-1)+time|id,
               method='REML',
               data=y_mi,
               control=lmeControl(returnObject=TRUE))
    ana <- summary(lmi)$coef$fixed
    d_ana <- c(0, bl_T-bl_P, 0, -l/N, 1-r/N, r/N, -l/N, 1-r/N, r/N, bl_T-bl_P)
    theta_mar_mi_j[j] <- d_ana%*%ana
    vartotal_mar_mi_j[j] <- d_ana%*%vcov(lmi)%*%d_ana
    
    # Mild
    y_mis$t2[ind_mis] <- NA
    y_misp$t2[indp_mis] <- NA
    for(k in 1:N){
      if(is.na(y_mis$t2[k])){
        z_act <- matrix(c(1,y_mis$t0[k],0,0,1,0,0,0,0,0,
                          1,y_mis$t0[k],1,0,0,0,0,0,0,y_mis$t0[k]),nrow=2,byrow=T)
        mu_act <- z_act %*% coeff
        mu_act_tilde <- mvrnorm(1, mu_act, sigma)
        y_mis$t2[k] <- rnorm(1,
                             mu_act_tilde[2]+sigma_tilde[2,1]/sigma_tilde[1,1]*(y_mis$t1[k]-mu_act_tilde[1]),
                             sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
      }
      if(is.na(y_misp$t2[k])){
        z_act <- matrix(c(1,y_misp$t0[k],0,0,0,0,0,0,0,0,
                          1,y_misp$t0[k],1,0,0,0,0,0,0,y_misp$t0[k]),nrow=2,byrow=T)
        mu_act <- z_act %*% coeff
        mu_act_tilde <- mvrnorm(1, mu_act, sigma)
        y_misp$t2[k] <- rnorm(1,
                              mu_act_tilde[2]+sigma_tilde[2,1]/sigma_tilde[1,1]*(y_misp$t1[k]-mu_act_tilde[1]),
                              sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
      }
    }
    y_mi_tot <- rbind(y_mis,y_misp)
    y_mi <- melt(y_mi_tot,id.vars = c('id','t0','arm'))
    colnames(y_mi) <- c('id','bl','arm','time','val')
    lmi <- lme(val~bl+time+arm+arm:time+bl:time,
               random=~(-1)+time|id,
               method='REML',
               data=y_mi,
               control=lmeControl(returnObject=TRUE))
    ana <- summary(lmi)$coef$fixed
    d_ana <- c(0, bl_T-bl_P, 0, -l/N, 1-r/N, r/N, -l/N, 1-r/N, r/N, bl_T-bl_P)
    theta_mild_mi_j[j] <- d_ana%*%ana
    vartotal_mild_mi_j[j] <- d_ana%*%vcov(lmi)%*%d_ana
    
    # Moderate
    y_mis$t2[ind_mis] <- NA
    y_misp$t2[indp_mis] <- NA
    for(k in 1:N){
      if(is.na(y_mis$t2[k])){
        z_act <- matrix(c(1,y_mis$t0[k],0,0,1,0,0,0,0,0,
                          1,y_mis$t0[k],1,0,0,0,0,0,0,y_mis$t0[k]),nrow=2,byrow=T)
        mu_act <- z_act %*% coeff
        mu_act_tilde <- mvrnorm(1, mu_act, sigma)
        y_mis$t2[k] <- rnorm(1,
                             y_mis$t0[k]+sigma_tilde[2,1]/sigma_tilde[1,1]*(y_mis$t1[k]-mu_act_tilde[1]),
                             sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
      }
      if(is.na(y_misp$t2[k])){
        z_act <- matrix(c(1,y_misp$t0[k],0,0,0,0,0,0,0,0,
                          1,y_misp$t0[k],1,0,0,0,0,0,0,y_misp$t0[k]),nrow=2,byrow=T)
        mu_act <- z_act %*% coeff
        mu_act_tilde <- mvrnorm(1, mu_act, sigma)
        y_misp$t2[k] <- rnorm(1,
                              mu_act_tilde[2]+sigma_tilde[2,1]/sigma_tilde[1,1]*(y_misp$t1[k]-mu_act_tilde[1]),
                              sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
      }
    }
    y_mi_tot <- rbind(y_mis,y_misp)
    y_mi <- melt(y_mi_tot,id.vars = c('id','t0','arm'))
    colnames(y_mi) <- c('id','bl','arm','time','val')
    lmi <- lme(val~bl+time+arm+arm:time+bl:time,
               random=~(-1)+time|id,
               method='REML',
               data=y_mi,
               control=lmeControl(returnObject=TRUE))
    ana <- summary(lmi)$coef$fixed
    theta_moderate_mi_j[j] <- d_ana%*%ana
    vartotal_moderate_mi_j[j] <- d_ana%*%vcov(lmi)%*%d_ana
    
    # Extreme
    y_mis$t2[ind_mis] <- NA
    y_misp$t2[indp_mis] <- NA
    for(k in 1:N){
      if(is.na(y_mis$t2[k])){
        z_act <- matrix(c(1,y_mis$t0[k],0,0,1,0,0,0,0,0,
                          1,y_mis$t0[k],1,0,1,0,0,1,0,y_mis$t0[k]),nrow=2,byrow=T)
        mu_act <- z_act %*% coeff
        mu_act_tilde <- mvrnorm(1, mu_act, sigma)
        y_mis$t2[k] <- rnorm(1,
                             y_mis$t0[k]+sigma_tilde[2,1]/sigma_tilde[1,1]*(y_mis$t1[k]-mu_act_tilde[1]),
                             sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
      }
      if(is.na(y_misp$t2[k])){
        z_act <- matrix(c(1,y_misp$t0[k],0,0,0,0,0,0,0,0,
                          1,y_misp$t0[k],1,0,0,1,0,0,1,y_misp$t0[k]),nrow=2,byrow=T)
        mu_act <- z_act %*% coeff
        mu_act_tilde <- mvrnorm(1, mu_act, sigma)
        y_misp$t2[k] <- rnorm(1,
                              mu_act_tilde[2]+sigma_tilde[2,1]/sigma_tilde[1,1]*(y_misp$t1[k]-mu_act_tilde[1]),
                              sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
        
      }
    }
    y_mi_tot <- rbind(y_mis,y_misp)
    y_mi <- melt(y_mi_tot,id.vars = c('id','t0','arm'))
    colnames(y_mi) <- c('id','bl','arm','time','val')
    lmi <- lme(val~bl+time+arm+arm:time+bl:time,
               random=~(-1)+time|id,
               method='REML',
               data=y_mi,
               control=lmeControl(returnObject=TRUE))
    ana <- summary(lmi)$coef$fixed
    theta_extreme_mi_j[j] <- d_ana%*%ana
    vartotal_extreme_mi_j[j] <- d_ana%*%vcov(lmi)%*%d_ana
  }
  # MAR
  theta_mar_mi <- mean(theta_mar_mi_j)
  W <- mean(vartotal_mar_mi_j)
  B <- var(theta_mar_mi_j)
  se_mar_mi <-  sqrt(W + (1+1/MI)*B)
  df <- (1+MI/(MI+1)*W/B)*(MI-1)
  low_mar_mi <- theta_mar_mi - qt(1-alpha/2, df)*se_mar_mi
  upp_mar_mi <- theta_mar_mi + qt(1-alpha/2, df)*se_mar_mi
  cover_mar_mi <- low_mar_mi<true_df & true_df<upp_mar_mi
  t <- theta_mar_mi/se_mar_mi
  pval_mar_mi <- (1-pt(t, df)<.025)*100
  out <- c(out, theta_mar_mi, se_mar_mi, cover_mar_mi, pval_mar_mi)
  # Mild
  theta_mild_mi <- mean(theta_mild_mi_j)
  W <- mean(vartotal_mild_mi_j)
  B <- var(theta_mild_mi_j)
  se_mild_mi <-  sqrt(W + (1+1/MI)*B)
  df <- (1+MI/(MI+1)*W/B)*(MI-1)
  low_mild_mi <- theta_mild_mi - qt(1-alpha/2, df)*se_mild_mi
  upp_mild_mi <- theta_mild_mi + qt(1-alpha/2, df)*se_mild_mi
  cover_mild_mi <- low_mild_mi<true_df & true_df<upp_mild_mi
  t <- theta_mild_mi/se_mild_mi
  pval_mild_mi <- (1-pt(t, df)<.025)*100
  out <- c(out, theta_mild_mi, se_mild_mi, cover_mild_mi, pval_mild_mi)
  # Moderate
  theta_moderate_mi <- mean(theta_moderate_mi_j)
  W <- mean(vartotal_moderate_mi_j)
  B <- var(theta_moderate_mi_j)
  se_moderate_mi <-  sqrt(W + (1+1/MI)*B)
  df <- (1+MI/(MI+1)*W/B)*(MI-1)
  low_moderate_mi <- theta_moderate_mi - qt(1-alpha/2, df)*se_moderate_mi
  upp_moderate_mi <- theta_moderate_mi + qt(1-alpha/2, df)*se_moderate_mi
  cover_moderate_mi <- low_moderate_mi<true_df & true_df<upp_moderate_mi
  t <- theta_moderate_mi/se_moderate_mi
  pval_moderate_mi <- (1-pt(t, df)<.025)*100
  out <- c(out, theta_moderate_mi, se_moderate_mi, cover_moderate_mi, pval_moderate_mi)
  # Extreme
  theta_extreme_mi <- mean(theta_extreme_mi_j)
  W <- mean(vartotal_extreme_mi_j)
  B <- var(theta_extreme_mi_j)
  se_extreme_mi <-  sqrt(W + (1+1/MI)*B)
  df <- (1+MI/(MI+1)*W/B)*(MI-1)
  low_extreme_mi <- theta_extreme_mi - qt(1-alpha/2, df)*se_extreme_mi
  upp_extreme_mi <- theta_extreme_mi + qt(1-alpha/2, df)*se_extreme_mi
  cover_extreme_mi <- low_extreme_mi<true_df & true_df<upp_extreme_mi
  t <- theta_extreme_mi/se_extreme_mi
  pval_extreme_mi <- (1-pt(t, df)<.025)*100
  out <- c(out, theta_extreme_mi, se_extreme_mi, cover_extreme_mi, pval_extreme_mi)
  
  # Bayesian Regression
  MI <- 20
  theta_mar_reg_j <- vartotal_mar_reg_j <-
    theta_mild_reg_j <- vartotal_mild_reg_j <-
    theta_moderate_reg_j <- vartotal_moderate_reg_j <-
    theta_extreme_reg_j <- vartotal_extreme_reg_j <-
    rep(NA, MI)
  y_mis$t2[ind_mis] <- NA
  y_misp$t2[indp_mis] <- NA
  y_reg <- rbind(y_mis,y_misp)
  y_reg$trt <- c(rep(1,N),rep(0,N))
  y_reg$trt[ind_df] <- 0
  y_reg$tact <- rep(0,2*N)
  y_reg$tact[ind_df] <- 1
  y_reg$pla <- c(rep(0,N),rep(1,N))
  y_reg$pla[N+indp_df] <- 0
  y_reg$pact <- rep(0,2*N)
  y_reg$pact[N+indp_df] <- 1
  
  imputed <- mice(y_reg, m=MI, method='norm')$imp$t2
  
  bl_DT <- mean(y_mis$t0[-ind_df])
  bl_RT <- mean(y_mis$t0[ind_res])
  bl_MT <- mean(y_mis$t0[ind_mis])
  bl_DP <- mean(y_misp$t0[-ind_df])
  bl_RP <- mean(y_misp$t0[ind_res])
  bl_MP <- mean(y_misp$t0[ind_mis])
  
  dt_mild <- coeff[5]+coeff[8]
  dt_moderate <- coeff[1]+bl_MT*coeff[2]+coeff[3]+coeff[5]+coeff[8]+bl_MT*coeff[10]-bl_MT
  dt_extreme <- coeff[6]+coeff[9]
  y_reg <- y_reg[,1:5]
  for(j in 1:MI){
    # MAR
    y_reg$t2[c(ind_mis, N+indp_mis)] <- imputed[,j]
    y_mi <- melt(y_reg,id.vars = c('id','t0','arm'))
    colnames(y_mi) <- c('id','bl','arm','time','val')
    lmi <- lme(val~bl+time+arm+arm:time+bl:time,
               random=~(-1)+time|id,
               method='REML',
               data=y_mi,
               control=lmeControl(returnObject=TRUE))
    ana <- summary(lmi)$coef$fixed
    d_ana <- c(0, bl_T-bl_P, 0, -l/N, 1-r/N, r/N, -l/N, 1-r/N, r/N, bl_T-bl_P)
    theta_mar_reg_j[j] <- d_ana%*%ana
    vartotal_mar_reg_j[j] <- d_ana%*%vcov(lmi)%*%d_ana
    # Mild
    y_reg$t2[c(ind_mis, N+indp_mis)] <- imputed[,j]-c(rep(dt_mild, length(ind_mis)),rep(0, length(indp_mis)))
    y_mi <- melt(y_reg,id.vars = c('id','t0','arm'))
    colnames(y_mi) <- c('id','bl','arm','time','val')
    lmi <- lme(val~bl+time+arm+arm:time+bl:time,
               random=~(-1)+time|id,
               method='REML',
               data=y_mi,
               control=lmeControl(returnObject=TRUE))
    ana <- summary(lmi)$coef$fixed
    d_ana <- c(0, bl_T-bl_P, 0, -l/N, 1-r/N, r/N, -l/N, 1-r/N, r/N, bl_T-bl_P)
    
    theta_mild_reg_j[j] <- d_ana%*%ana
    vartotal_mild_reg_j[j] <- d_ana%*%vcov(lmi)%*%d_ana
    # Moderate
    y_reg$t2[c(ind_mis, N+indp_mis)] <- imputed[,j]-c(rep(dt_moderate, length(ind_mis)),rep(0, length(indp_mis)))
    y_mi <- melt(y_reg,id.vars = c('id','t0','arm'))
    colnames(y_mi) <- c('id','bl','arm','time','val')
    lmi <- lme(val~bl+time+arm+arm:time+bl:time,
               random=~(-1)+time|id,
               method='REML',
               data=y_mi,
               control=lmeControl(returnObject=TRUE))
    ana <- summary(lmi)$coef$fixed
    d_ana <- c(0, bl_T-bl_P, 0, -l/N, 1-r/N, r/N, -l/N, 1-r/N, r/N, bl_T-bl_P)
    theta_moderate_reg_j[j] <- d_ana%*%ana
    vartotal_moderate_reg_j[j] <- d_ana%*%vcov(lmi)%*%d_ana
    # Extreme
    y_reg$t2[c(ind_mis, N+indp_mis)] <- imputed[,j]-c(rep(dt_moderate, length(ind_mis)),rep(-dt_extreme, length(indp_mis)))
    y_mi <- melt(y_reg,id.vars = c('id','t0','arm'))
    colnames(y_mi) <- c('id','bl','arm','time','val')
    lmi <- lme(val~bl+time+arm+arm:time+bl:time,
               random=~(-1)+time|id,
               method='REML',
               data=y_mi,
               control=lmeControl(returnObject=TRUE))
    ana <- summary(lmi)$coef$fixed
    d_ana <- c(0, bl_T-bl_P, 0, -l/N, 1-r/N, r/N, -l/N, 1-r/N, r/N, bl_T-bl_P)
    theta_extreme_reg_j[j] <- d_ana%*%ana
    vartotal_extreme_reg_j[j] <- d_ana%*%vcov(lmi)%*%d_ana
  }
  # MAR
  theta_mar_reg <- mean(theta_mar_reg_j)
  W <- mean(vartotal_mar_reg_j)
  B <- var(theta_mar_reg_j)
  se_mar_reg <-  sqrt(W + (1+1/MI)*B)
  df <- (1+MI/(MI+1)*W/B)*(MI-1)
  low_mar_reg <- theta_mar_reg - qt(1-alpha/2, df)*se_mar_reg
  upp_mar_reg <- theta_mar_reg + qt(1-alpha/2, df)*se_mar_reg
  cover_mar_reg <- low_mar_reg<true_df & true_df<upp_mar_reg
  t <- theta_mar_reg/se_mar_reg
  pval_mar_reg <- (1-pt(t, df)<.025)*100
  out <- c(out, theta_mar_reg, se_mar_reg, cover_mar_reg, pval_mar_reg)
  # Mild
  theta_mild_reg <- mean(theta_mild_reg_j)
  W <- mean(vartotal_mild_reg_j)
  B <- var(theta_mild_reg_j)
  se_mild_reg <-  sqrt(W + (1+1/MI)*B)
  df <- (1+MI/(MI+1)*W/B)*(MI-1)
  low_mild_reg <- theta_mild_reg - qt(1-alpha/2, df)*se_mild_reg
  upp_mild_reg <- theta_mild_reg + qt(1-alpha/2, df)*se_mild_reg
  cover_mild_reg <- low_mild_reg<true_df & true_df<upp_mild_reg
  t <- theta_mild_reg/se_mild_reg
  pval_mild_reg <- (1-pt(t, df)<.025)*100
  out <- c(out, theta_mild_reg, se_mild_reg, cover_mild_reg, pval_mild_reg)
  # Moderate
  theta_moderate_reg <- mean(theta_moderate_reg_j)
  W <- mean(vartotal_moderate_reg_j)
  B <- var(theta_moderate_reg_j)
  se_moderate_reg <-  sqrt(W + (1+1/MI)*B)
  df <- (1+MI/(MI+1)*W/B)*(MI-1)
  low_moderate_reg <- theta_moderate_reg - qt(1-alpha/2, df)*se_moderate_reg
  upp_moderate_reg <- theta_moderate_reg + qt(1-alpha/2, df)*se_moderate_reg
  cover_moderate_reg <- low_moderate_reg<true_df & true_df<upp_moderate_reg
  t <- theta_moderate_reg/se_moderate_reg
  pval_moderate_reg <- (1-pt(t, df)<.025)*100
  out <- c(out, theta_moderate_reg, se_moderate_reg, cover_moderate_reg, pval_moderate_reg)
  # Extreme
  theta_extreme_reg <- mean(theta_extreme_reg_j)
  W <- mean(vartotal_extreme_reg_j)
  B <- var(theta_extreme_reg_j)
  se_extreme_reg <-  sqrt(W + (1+1/MI)*B)
  df <- (1+MI/(MI+1)*W/B)*(MI-1)
  low_extreme_reg <- theta_extreme_reg - qt(1-alpha/2, df)*se_extreme_reg
  upp_extreme_reg <- theta_extreme_reg + qt(1-alpha/2, df)*se_extreme_reg
  cover_extreme_reg <- low_extreme_reg<true_df & true_df<upp_extreme_reg
  t <- theta_extreme_reg/se_extreme_reg
  pval_extreme_reg <- (1-pt(t, df)<.025)*100
  out <- c(out, theta_extreme_reg, se_extreme_reg, cover_extreme_reg, pval_extreme_reg)
})

stopCluster(cl)

round(rowMeans(result),2)
round(sqrt(apply(result,1,var)),2)
