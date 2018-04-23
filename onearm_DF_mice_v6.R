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

S <- 100

df <- switch <- miss <- c(NA,S)
for(i in 1:S){
  y_dj <- mvrnorm(N, mu_dj, vv)
  y_dj <- data.frame(id=1:N, t0=y_dj[,1], t1=y_dj[,2], t2=y_dj[,3], arm='trt',stringsAsFactors=FALSE)
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
  df[i] <- mean(y_df$t2)
  px <- as.matrix(cbind(1,y_df[ind_df,3:4])) %*% t(t(psi2))
  pt <- 1/(1+exp(-px))
  ind[ind==2] <- ind[ind==2] + rbinom(sum(ind==2), 1, pt)
  while (sum(ind==3)==0|sum(ind==2)==0){
    ind[ind==3] <- 2
    ind[ind==2]<-ind[ind==2] + rbinom(sum(ind==2), 1, pt)
  } 
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
  library(emmeans)
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
  while (sum(ind==3)==0|sum(ind==2)==0) {
    ind[ind==3] <- 2
    ind[ind==2]<-ind[ind==2] + rbinom(sum(ind==2), 1, pt)
  }
  ind_mis <- which(ind==3)
  ind_res <- which(ind==2)
  y_mis <- y_df
  y_mis$t2[ind_mis] <- NA
  y_mis_long <- melt(y_mis,id.vars = c('id','trt','res','t0'))
  colnames(y_mis_long) <- c('id','trt','res','bl','time','val')
  y_mis_long <- na.omit(y_mis_long)
  
  r <- N-length(ind_df)
  q <- N-r-length(ind_mis)
  
  
  # MMRM
  m1 <- lme(val~bl+res+time+res:time+bl:time,
            random=~(-1)+time|id,
            method='REML',
            data=y_mis_long,
            control=lmeControl(returnObject=TRUE))
  bl_bar <- mean(y_df$t0)
  bl_O <- mean(y_df$t0[-ind_mis])
  bl_D <- mean(y_df$t0[ind_df])
  bl_R <- mean(y_df$t0[ind_res])
  bl_M <- mean(y_df$t0[ind_mis])
  # MAR-MMRM
  coeff <- summary(m1)$coef$fixed
  D_mar <- c(1, bl_bar, 1-r/N, 1, 1-r/N, bl_bar)
  theta_mar_mmrm <- D_mar%*%coeff
  se_mar_mmrm <- sqrt(D_mar%*%vcov(m1)%*%D_mar + (coeff[3]+coeff[5])^2*r*(N-r)/N^3)
  low_mar_mmrm <- theta_mar_mmrm - qt(1-alpha/2, N-2)*se_mar_mmrm
  upp_mar_mmrm <- theta_mar_mmrm + qt(1-alpha/2, N-2)*se_mar_mmrm
  cover_mar_mmrm <- low_mar_mmrm<true_df & true_df<upp_mar_mmrm
  out <- c(out, theta_mar_mmrm, se_mar_mmrm, cover_mar_mmrm)
  # MTN-MMRM
  D_mtn <- c(1, bl_bar, (N-r)/N, (r+q)/N, q/N, (r+q)/N*bl_O)
  theta_mtn_mmrm <- D_mtn%*%coeff
  se_mtn_mmrm <- sqrt(D_mtn%*%vcov(m1)%*%D_mtn + 
                        r*(N-r)/N^3*(coeff[1]+(coeff[2]+coeff[6])*bl_D+coeff[4])^2+
                        q*(N-q)/N^3*(coeff[1]+(coeff[2]+coeff[6])*bl_R+coeff[3]+coeff[4]+coeff[5])^2+
                        (r+q)*(N-r-q)/N^3*(coeff[1]+coeff[2]*bl_M+coeff[3])^2-
                        2*r*q/N^3*(coeff[1]+(coeff[2]+coeff[6])*bl_D+coeff[4])*(coeff[1]+(coeff[2]+coeff[6])*bl_R+coeff[3]+coeff[4]+coeff[5])-
                        2*r*(N-r-q)/N^3*(coeff[1]+(coeff[2]+coeff[6])*bl_D+coeff[4])*(coeff[1]+coeff[2]*bl_M+coeff[3])-
                        2*q*(N-r-q)/N^3*(coeff[1]+(coeff[2]+coeff[6])*bl_R+coeff[3]+coeff[4]+coeff[5])*(coeff[1]+coeff[2]*bl_M+coeff[3])
  )
  low_mtn_mmrm <- theta_mtn_mmrm - qt(1-alpha/2, N-2)*se_mtn_mmrm
  upp_mtn_mmrm <- theta_mtn_mmrm + qt(1-alpha/2, N-2)*se_mtn_mmrm
  cover_mtn_mmrm <- low_mtn_mmrm<true_df & true_df<upp_mtn_mmrm
  out <- c(out, theta_mtn_mmrm, se_mtn_mmrm, cover_mtn_mmrm)
  # J2B-MMRM
  D_j2b <- c((r+q)/N, (r+q)*bl_O/N, q/N, (r+q)/N, q/N, (r+q)*bl_O/N)
  theta_j2b_mmrm <- D_j2b%*%coeff+(N-r-q)/N*bl_M
  se_j2b_mmrm <- sqrt(D_j2b%*%vcov(m1)%*%D_j2b +
                        (coeff[1]+bl_D*(coeff[2]+coeff[6])+coeff[4])^2*r*(N-r)/N^3+
                        (coeff[1]+bl_R*(coeff[2]+coeff[6])+coeff[3]+coeff[4]+coeff[5])^2*q*(N-q)/N^3+
                        bl_M^2*(N-r-q)*(r+q)/N^3-
                        2*r*q/N^3*(coeff[1]+bl_D*(coeff[2]+coeff[6])+coeff[4])*(coeff[1]+bl_R*(coeff[2]+coeff[6])+coeff[3]+coeff[4]+coeff[5])-
                        2*r*(N-r-q)/N^3*(coeff[1]+bl_D*(coeff[2]+coeff[6])+coeff[4])*bl_M-
                        2*q*(N-r-q)/N^3*(coeff[1]+bl_R*(coeff[2]+coeff[6])+coeff[3]+coeff[4]+coeff[5])*bl_M
  )
  low_j2b_mmrm <- theta_j2b_mmrm - qt(1-alpha/2, N-2)*se_j2b_mmrm
  upp_j2b_mmrm <- theta_j2b_mmrm + qt(1-alpha/2, N-2)*se_j2b_mmrm
  cover_j2b_mmrm <- low_j2b_mmrm<true_df & true_df<upp_j2b_mmrm
  out <- c(out, theta_j2b_mmrm, se_j2b_mmrm, cover_j2b_mmrm)
  # MI
  MI <- 20
  sigma <- getVarCov(m1) + diag(2)*m1$sigma^2
  theta_mar_mi_j <- vartotal_mar_mi_j <- 
    theta_mtn_mi_j <- vartotal_mtn_mi_j <-
    theta_j2b_mi_j <- vartotal_j2b_mi_j <-
    rep(NA, MI)
  for(j in 1:MI){
    # MAR
    y_mis$t2[ind_mis] <- NA
    for(k in 1:N){
      if(is.na(y_mis$t2[k])){
        z_act <- matrix(c(1,y_mis$t0[k],1,0,0,0,
                          1,y_mis$t0[k],1,1,1,y_mis$t0[k]),nrow=2,byrow=T)
        mu_act <- z_act %*% coeff
        mu_act_tilde <- mvrnorm(1, mu_act, sigma)
        sigma_tilde <- matrix(rWishart(1, N-2, sigma),nrow=2)/(N-2)
        y_mis$t2[k] <- rnorm(1,
                             mu_act_tilde[2]+sigma_tilde[2,1]/sigma_tilde[1,1]*(y_mis$t1[k]-mu_act_tilde[1]),
                             sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
      }
    }
    y_mi <- melt(y_mis,id.vars = c('id','trt','res','t0'))
    colnames(y_mi) <- c('id','trt','res','bl','time','val')
    l <- summary(lsmeans(lme(val~bl*time,
                             random=~(-1)+time|id,
                             method='REML',
                             data=y_mi,
                             control=lmeControl(returnObject=TRUE)), ~time, data=y_mi))[2,]
    theta_mar_mi_j[j] <- l$lsmean
    vartotal_mar_mi_j[j] <- l$SE^2
    
    # MTN
    y_mis$t2[ind_mis] <- NA
    for(k in 1:N){
      if(is.na(y_mis$t2[k])){
        z_act <- matrix(c(1,y_mis$t0[k],1,0,0,0,
                          1,y_mis$t0[k],1,1,1,y_mis$t0[k]),nrow=2,byrow=T)
        mu_act <- z_act %*% coeff
        mu_act_tilde <- mvrnorm(1, mu_act, sigma)
        sigma_tilde <- matrix(rWishart(1, N-2, sigma),nrow=2)/(N-2)
        y_mis$t2[k] <- rnorm(1,
                             mu_act_tilde[1]+sigma_tilde[2,1]/sigma_tilde[1,1]*(y_mis$t1[k]-mu_act_tilde[1]),
                             sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
      }
    }
    y_mi <- melt(y_mis,id.vars = c('id','trt','res','t0'))
    colnames(y_mi) <- c('id','trt','res','bl','time','val')
    l <- summary(lsmeans(lme(val~bl*time,
                             random=~(-1)+time|id,
                             method='REML',
                             data=y_mi,
                             control=lmeControl(returnObject=TRUE)), ~time, data=y_mi))[2,]
    theta_mtn_mi_j[j] <- l$lsmean
    vartotal_mtn_mi_j[j] <- l$SE^2
    # J2B
    y_mis$t2[ind_mis] <- NA
    for(k in 1:N){
      if(is.na(y_mis$t2[k])){
        z_act <- matrix(c(1,y_mis$t0[k],1,0,0,0,
                          1,y_mis$t0[k],1,1,1,y_mis$t0[k]),nrow=2,byrow=T)
        mu_act <- z_act %*% coeff
        mu_act_tilde <- mvrnorm(1, mu_act, sigma)
        sigma_tilde <- matrix(rWishart(1, N-2, sigma),nrow=2)/(N-2)
        y_mis$t2[k] <- rnorm(1,
                             y_mis$t0[k]+sigma_tilde[2,1]/sigma_tilde[1,1]*(y_mis$t1[k]-mu_act_tilde[1]),
                             sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
        
      }
    }
    y_mi <- melt(y_mis,id.vars = c('id','trt','res','t0'))
    colnames(y_mi) <- c('id','trt','res','bl','time','val')
    l <- summary(lsmeans(lme(val~bl*time,
                             random=~(-1)+time|id,
                             method='REML',
                             data=y_mi,
                             control=lmeControl(returnObject=TRUE)), ~time, data=y_mi))[2,]
    theta_j2b_mi_j[j] <- l$lsmean
    vartotal_j2b_mi_j[j] <- l$SE^2
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
  out <- c(out, theta_mar_mi, se_mar_mi, cover_mar_mi)
  # MTN
  theta_mtn_mi <- mean(theta_mtn_mi_j)
  W <- mean(vartotal_mtn_mi_j)
  B <- var(theta_mtn_mi_j)
  se_mtn_mi <-  sqrt(W + (1+1/MI)*B)
  df <- (1+MI/(MI+1)*W/B)*(MI-1)
  low_mtn_mi <- theta_mtn_mi - qt(1-alpha/2, df)*se_mtn_mi
  upp_mtn_mi <- theta_mtn_mi + qt(1-alpha/2, df)*se_mtn_mi
  cover_mtn_mi <- low_mtn_mi<true_df & true_df<upp_mtn_mi
  out <- c(out, theta_mtn_mi, se_mtn_mi, cover_mtn_mi)
  # J2B
  theta_j2b_mi <- mean(theta_j2b_mi_j)
  W <- mean(vartotal_j2b_mi_j)
  B <- var(theta_j2b_mi_j)
  se_j2b_mi <-  sqrt(W + (1+1/MI)*B)
  df <- (1+MI/(MI+1)*W/B)*(MI-1)
  low_j2b_mi <- theta_j2b_mi - qt(1-alpha/2, df)*se_j2b_mi
  upp_j2b_mi <- theta_j2b_mi + qt(1-alpha/2, df)*se_j2b_mi
  cover_j2b_mi <- low_j2b_mi<true_df & true_df<upp_j2b_mi
  out <- c(out, theta_j2b_mi, se_j2b_mi, cover_j2b_mi)
  
  # Sequential Regression
  MI <- 20
  theta_mar_reg_j <- vartotal_mar_reg_j <-
    theta_mtn_reg_j <- vartotal_mtn_reg_j <-
    theta_j2b_reg_j <- vartotal_j2b_reg_j <-
    rep(NA, MI)
  y_mis$t2[ind_mis] <- NA
  imputed <- mice(y_mis, m=MI, method='norm')$imp$t2
  l <- summary(lm(t2~t0+t1+res,data=y_mis))$coeff[,1]
  l2 <- summary(lm(t1~t0+res,data=y_mis))$coeff[,1]
  dt_mtn <- l[1]+l[2]*bl_M+l[3]*mean(y_df$t1[ind_mis])+l[4]-l2[1]-l2[2]*bl_M-l2[3]
  dt_j2b <- l[1]+(l[2]-1)*bl_M+l[3]*mean(y_df$t1[ind_mis])+l[4]
  
  for(j in 1:MI){
    # MAR
    y_mar <- y_mis
    y_mar$t2[ind_mis] <- imputed[,j]
    y_mar_long <- melt(y_mar,id.vars = c('id','trt','res','t0'))
    colnames(y_mar_long) <- c('id','trt','res','bl','time','val')
    lmi <- summary(lsmeans(lme(val~bl*time,
                               random=~(-1)+time|id,
                               method='REML',
                               data=y_mar_long,
                               control=lmeControl(returnObject=TRUE)), ~time, data=y_mar_long))[2,]
    theta_mar_reg_j[j] <- lmi$lsmean
    vartotal_mar_reg_j[j] <- lmi$SE^2
    
    # MTN
    y_mtn <- y_mis
    y_mtn$t2[ind_mis] <- imputed[,j]-dt_mtn
    y_mtn_long <- melt(y_mtn,id.vars = c('id','trt','res','t0'))
    colnames(y_mtn_long) <- c('id','trt','res','bl','time','val')
    lmi <- summary(lsmeans(lme(val~bl*time,
                               random=~(-1)+time|id,
                               method='REML',
                               data=y_mtn_long,
                               control=lmeControl(returnObject=TRUE)), ~time, y_mtn_long))[2,]
    theta_mtn_reg_j[j] <- lmi$lsmean
    vartotal_mtn_reg_j[j] <- lmi$SE^2
    
    # J2B
    y_j2b <- y_mis
    y_j2b$t2[ind_mis] <- imputed[,j]-dt_j2b
    y_j2b_long <- melt(y_j2b,id.vars = c('id','trt','res','t0'))
    colnames(y_j2b_long) <- c('id','trt','res','bl','time','val')
    lmi <- summary(lsmeans(lme(val~bl*time,
                               random=~(-1)+time|id,
                               method='REML',
                               data=y_j2b_long,
                               control=lmeControl(returnObject=TRUE)), ~time, y_j2b_long))[2,]
    theta_j2b_reg_j[j] <- lmi$lsmean
    vartotal_j2b_reg_j[j] <- lmi$SE^2
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
  out <- c(out, theta_mar_reg, se_mar_reg, cover_mar_reg)
  
  # MTN
  theta_mtn_reg <- mean(theta_mtn_reg_j)
  W <- mean(vartotal_mtn_reg_j)
  B <- var(theta_mtn_reg_j)
  se_mtn_reg <-  sqrt(W + (1+1/MI)*B)
  df <- (1+MI/(MI+1)*W/B)*(MI-1)
  low_mtn_reg <- theta_mtn_reg - qt(1-alpha/2, df)*se_mtn_reg
  upp_mtn_reg <- theta_mtn_reg + qt(1-alpha/2, df)*se_mtn_reg
  cover_mtn_reg <- low_mtn_reg<true_df & true_df<upp_mtn_reg
  out <- c(out, theta_mtn_reg, se_mtn_reg, cover_mtn_reg)
  
  # J2B
  theta_j2b_reg <- mean(theta_j2b_reg_j)
  W <- mean(vartotal_j2b_reg_j)
  B <- var(theta_j2b_reg_j)
  se_j2b_reg <-  sqrt(W + (1+1/MI)*B)
  df <- (1+MI/(MI+1)*W/B)*(MI-1)
  low_j2b_reg <- theta_j2b_reg - qt(1-alpha/2, df)*se_j2b_reg
  upp_j2b_reg <- theta_j2b_reg + qt(1-alpha/2, df)*se_j2b_reg
  cover_j2b_reg <- low_j2b_reg<true_df & true_df<upp_j2b_reg
  out <- c(out, theta_j2b_reg, se_j2b_reg, cover_j2b_reg)
})

stopCluster(cl)

round(rowMeans(result),2)
round(sqrt(apply(result,1,var)),2)
