library(MASS)
library(nlme)
library(reshape2)
library(parallel)


rm(list = ls())


mu_dj <- c(0, 1, 4)
mu_df <- c(0, 1, 2)

sd <- c(2.0, 1.8, 2.2)
corr <- matrix(c(1.0, 0.6, 0.1,
                 0.6, 1.0, 0.2,
                 0.1, 0.2, 1.0), nrow=3)
vv <- diag(sd) %*% corr %*% diag(sd)
N <- 100
alpha <- .05

# S1
psi1 <- c(0.4, -0.4, 0)
psi2 <- c(-0.3, -0.4, 0)
# S2
# psi1 <- c(0.4, -0.4, 0)
# psi2 <- c(0.3, 0, -0.4)
# S3
# psi1 <- c(1.6, 0, -0.4)
# psi2 <- c(-0.1, -0.4, 0)
# S4
# psi1 <- c(1.6, 0, -0.4)
# psi2 <- c(0.3, 0, -0.4)

S <- 5000

df <- switch <- miss <- c(NA,S)
for(i in 1:S){
  y_dj <- mvrnorm(N, mu_dj, vv)
  y_dj <- data.frame(id=1:N, t0=y_dj[,1], t1=y_dj[,2], t2=y_dj[,3],stringsAsFactors=FALSE)
  px <- as.matrix(cbind(1,y_dj[,3:4])) %*% t(t(psi1))
  pt <- 1/(1+exp(-px))
  ind <- rbinom(N, 1, pt)+1
  if (sum(ind==2)==0) ind <- rbinom(N, 1, pt)+1
  y_df <- y_dj
  ind_df <- which(ind==2)
  y_df$t2[ind_df] <- rnorm(length(ind_df),
                           mean=mu_df[3]+vv[3,1:2]%*%solve(vv[1:2,1:2])%*%t(sweep(y_df[ind_df,2:3], 2, mu_dj[1:2])),
                           sd=sqrt(vv[3,3]-vv[3,1:2]%*%solve(vv[1:2,1:2])%*%vv[1:2,3]))
  switch[i] <- length(ind_df)/N
  y_df$res <- 0
  y_df$res[ind_df] <- 1
  y_mis_long <- melt(y_df,id.vars = c('id','res','t0'))
  colnames(y_mis_long) <- c('id','res','bl','time','val')
  m1 <- lme(val~(-1)+bl+time,
            random=~(-1)+time|id,
            method='REML',
            data=y_mis_long,
            control=lmeControl(returnObject=TRUE))
  coeff <- summary(m1)$coef$fixed
  df[i] <- coeff[3]
  px <- as.matrix(cbind(1,y_df[ind_df,3:4])) %*% t(t(psi2))
  pt <- 1/(1+exp(-px))
  ind[ind==2] <- ind[ind==2] + rbinom(sum(ind==2), 1, pt)
  while (sum(ind==3)==0|sum(ind==2)==0) 
    ind[ind==2]<-ind[ind==2] + rbinom(sum(ind==2), 1, pt)
  ind_mis <- which(ind==3)
  miss[i] <- length(ind_mis)/length(ind_df)
}

true_df <- mean(df)
true_switch <- mean(switch)
true_miss <- mean(miss)
true_df
true_switch
true_miss




no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl,varlist=c('S', 'N', 'vv', 'alpha', 
                           'true_df',
                           'psi1', 'psi2', 
                           'mu_dj','mu_df'))
clusterSetRNGStream(cl, 1234)
clusterEvalQ(cl, {
  library(MASS)
  library(nlme)
  library(reshape2)
  library(mice)
  library(lsmeans)
})


result <- parSapply(cl, 1:S, function(x){
  out <- c()
  y_dj <- mvrnorm(N, mu_dj, vv)
  y_dj <- data.frame(id=1:N, t0=y_dj[,1], t1=y_dj[,2], t2=y_dj[,3],stringsAsFactors=FALSE)
  px <- as.matrix(cbind(1,y_dj[,3:4])) %*% t(t(psi1))
  pt <- 1/(1+exp(-px))
  ind <- rbinom(N, 1, pt)+1
  while (sum(ind==2)==0) ind <- rbinom(N, 1, pt)+1
  ind_df <- which(ind==2)
  y_df <- y_dj
  y_df$trt <- 1
  y_df$res <- 0
  y_df$trt[ind_df] <- 0
  y_df$res[ind_df] <- 1
  y_df$t2[ind_df] <- rnorm(length(ind_df), 
                           mean=mu_df[3]+vv[3,1:2]%*%solve(vv[1:2,1:2])%*%t(sweep(y_df[ind_df,2:3], 2, mu_dj[1:2])),
                           sd=sqrt(vv[3,3]-vv[3,1:2]%*%solve(vv[1:2,1:2])%*%vv[1:2,3]))
  px <- as.matrix(cbind(1,y_df[ind_df,3:4])) %*% t(t(psi2))
  pt <- 1/(1+exp(-px))
  ind[ind==2] <- ind[ind==2] + rbinom(sum(ind==2), 1, pt)
  while (sum(ind==3)==0|sum(ind==2)==0){
    ind[ind==3] <- 2
    ind[ind==2]<-ind[ind==2] + rbinom(sum(ind==2), 1, pt)
  } 
  ind_mis <- which(ind==3)
  y_mis <- y_df
  y_mis$t2[ind_mis] <- NA
  y_mis_long <- melt(y_mis,id.vars = c('id','trt','res','t0'))
  colnames(y_mis_long) <- c('id','trt','res','bl','time','val')
  y_mis_long <- na.omit(y_mis_long)
  
  r <- N-length(ind_df)
  q <- N-r-length(ind_mis)
  
  
  # MMRM
  m1 <- lme(val~(-1)+bl+res+time+res:time,
            random=~(-1)+time|id,
            method='REML',
            data=y_mis_long,
            control=lmeControl(returnObject=TRUE))

  # MAR-MMRM
  coeff <- summary(m1)$coef$fixed
  D_mar <- c(0,  1-r/N, 0, 1, 1-r/N)
  theta_mar_mmrm <- D_mar%*%coeff
  se_mar_mmrm <- sqrt(D_mar%*%vcov(m1)%*%D_mar + (coeff[2]+coeff[5])^2*r*(N-r)/N^3)
  low_mar_mmrm <- theta_mar_mmrm - qt(1-alpha/2, N-2)*se_mar_mmrm
  upp_mar_mmrm <- theta_mar_mmrm + qt(1-alpha/2, N-2)*se_mar_mmrm
  cover_mar_mmrm <- low_mar_mmrm<true_df & true_df<upp_mar_mmrm
  out <- c(out, theta_mar_mmrm, se_mar_mmrm, cover_mar_mmrm)
  
  # # MI
  # MI <- 20
  # sigma <- getVarCov(m1) + diag(2)*m1$sigma^2
  # theta_mar_mi_j <- vartotal_mar_mi_j <- 
  #   theta_mtn_mi_j <- vartotal_mtn_mi_j <-
  #   theta_j2b_mi_j <- vartotal_j2b_mi_j <-
  #   rep(NA, MI)
  # for(j in 1:MI){
  #   # MAR
  #   y_mis$t2[ind_mis] <- NA
  #   for(k in 1:N){
  #     if(is.na(y_mis$t2[k])){
  #       z_act <- matrix(c(1,y_mis$t0[k],1,0,0,0,
  #                         1,y_mis$t0[k],1,1,1,y_mis$t0[k]),nrow=2,byrow=T)
  #       mu_act <- z_act %*% coeff
  #       mu_act_tilde <- mvrnorm(1, mu_act, sigma)
  #       sigma_tilde <- matrix(rWishart(1, N-2, sigma),nrow=2)/(N-2)
  #       y_mis$t2[k] <- rnorm(1,
  #                            mu_act_tilde[2]+sigma_tilde[2,1]/sigma_tilde[1,1]*(y_mis$t1[k]-mu_act_tilde[1]),
  #                            sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
  #     }
  #   }
  #   y_mi <- melt(y_mis,id.vars = c('id','trt','res','t0'))
  #   colnames(y_mi) <- c('id','trt','res','bl','time','val')
  #   l <- lme(val~bl+res+time+res:time+bl:time,
  #            random=~(-1)+time|id,
  #            method='REML',
  #            data=y_mi,
  #            control=lmeControl(returnObject=TRUE))
  #   ana <- summary(l)$coef$fixed
  #   d_ana <- c(1, bl_bar, 1-r/N, 1, 1-r/N, bl_bar)
  #   theta_mar_mi_j[j] <- d_ana%*%ana
  #   vartotal_mar_mi_j[j] <- d_ana%*%vcov(l)%*%d_ana
  # }
  # # MAR
  # theta_mar_mi <- mean(theta_mar_mi_j)
  # W <- mean(vartotal_mar_mi_j)
  # B <- var(theta_mar_mi_j)
  # se_mar_mi <-  sqrt(W + (1+1/MI)*B)
  # df <- (1+MI/(MI+1)*W/B)*(MI-1)
  # low_mar_mi <- theta_mar_mi - qt(1-alpha/2, df)*se_mar_mi
  # upp_mar_mi <- theta_mar_mi + qt(1-alpha/2, df)*se_mar_mi
  # cover_mar_mi <- low_mar_mi<true_df & true_df<upp_mar_mi
  # out <- c(out, theta_mar_mi, se_mar_mi, cover_mar_mi)
  
})

stopCluster(cl)

round(rowMeans(result),2)
round(sqrt(apply(result,1,var)),2)
