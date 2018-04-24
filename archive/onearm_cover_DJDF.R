library(MASS)
library(nlme)
library(reshape2)
library(parallel)
library(invgamma)
library(MCMCpack)
library(Hmisc)
library(gmodels)

rm(list = ls())


mu_dj <- c(1, 4)
mu_df <- c(1, 2)

vv <- matrix(c(4, 3, 3, 9),nrow=2)
N <- 100
alpha <- .05

# psi1 <- c(1, -0.1, -0.1)
# psi2 <- c(2, -0.1, -0.1)
psi1 <- c(-1, -0.1, 0)
psi2 <- c(-2, -0.1, 0)
# psi1 <- c(1, -0.1, 0)
# psi2 <- c(.5, -0.1, 0)
# px1 <- psi1 %*% c(1,mu_dj)
# px2 <- psi2 %*% c(1,mu_dj)
# p1 <- 1/(1+exp(-px1))
# p2 <- 1/(1+exp(-px2))
S <- 5000

df <- c(NA,S)
for(i in 1:S){
  y_dj <- mvrnorm(N, mu_dj, vv)
  y_dj <- data.frame(id=1:N, t1=y_dj[,1], t2=y_dj[,2], arm='trt',stringsAsFactors=FALSE)
  px <- data.frame(p1=as.matrix(cbind(1,y_dj[,2:3])) %*% t(t(psi1)),
                   p2=as.matrix(cbind(1,y_dj[,2:3])) %*% t(t(psi2)))
  pi <- 1/(1+exp(-px))
  pi[2] <- pi[1]*pi[2]
  pi[1] <- pi[1]-pi[2]
  pi <- data.frame(p0=1-pi$p1-pi$p2, pi)
  
  ind <- apply(pi, 1, function(x){sample(1:3, 1, replace=T, prob=x)})
  y_df <- y_dj
  
  ind_df <- which(ind==2 | ind==3)
  y_df$arm[ind_df] <- 'act'
  y_df$t2[ind_df] <- rnorm(length(ind_df), 
                           mean=mu_df[2]+vv[2,1]/vv[2,2]*(y_df$t1[ind_df]-mu_df[1]),
                           sd=sqrt(vv[2,2]-vv[2,1]/vv[1,1]*vv[2,1]))
  df[i] <- mean(y_df$t2)
}

true_dj <- mu_dj[2]
true_df <- mean(df)



no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl,varlist=c('S', 'N', 'vv', 'alpha', 
                           'true_dj','true_df',
                           'psi1', 'psi2', 
                           'mu_dj', 'mu_df'))
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
  pi[2] <- pi[1]*pi[2]
  pi[1] <- pi[1]-pi[2]
  pi <- data.frame(p0=1-pi$p1-pi$p2, pi)
  
  ind <- apply(pi, 1, function(x){sample(1:3, 1, replace=T, prob=x)})
  while(length(unique(ind))<3|min(table(ind))<3)
    ind <- apply(pi, 1, function(x){sample(1:3, 1, replace=T, prob=x)})
  y_df <- y_dj

  ind_df <- which(ind==2 | ind==3)
  y_df$arm[ind_df] <- 'act'
  y_df$t2[ind_df] <- rnorm(length(ind_df), 
                     mean=mu_df[2]+vv[2,1]/vv[2,2]*(y_df$t1[ind_df]-mu_df[1]),
                     sd=sqrt(vv[2,2]-vv[2,1]/vv[1,1]*vv[2,1]))

  
  y_mis <- y_df
  ind_mis <- which(ind==3)
  y_mis$arm[ind_mis] <- 'mis'
  y_mis$t2[ind_mis] <- NA
  
  
  r <- length(ind_df)
  q <- length(ind_mis)
  
  ## De Jure complete
  theta_dj <- mean(y_dj$t2)
  vartotal_dj <- var(y_dj$t2)
  se_dj <- sqrt(vartotal_dj/N)
  low_dj <- theta_dj - qt(1-alpha/2, N-1)*se_dj
  upp_dj <- theta_dj + qt(1-alpha/2, N-1)*se_dj
  cover_dj <- low_dj<true_dj & true_dj<upp_dj
  
  out <- c(out, theta_dj, se_dj, cover_dj)
  ## De Facto complete
  theta_df <- mean(y_df$t2)
  vartotal_df <- var(y_df$t2)
  se_df <- sqrt(vartotal_df/N)
  low_df <- theta_df - qt(1-alpha/2, N-1)*se_df
  upp_df <- theta_df + qt(1-alpha/2, N-1)*se_df
  cover_df <- low_df<true_df & true_df<upp_df
  out <- c(out, theta_df, se_df, cover_df)
  
  
  
  ## De Jure missing
  # MMRM
  y_mis_long <- melt(y_mis,id.vars = c('id','arm'))
  colnames(y_mis_long) <- c('id','arm','time','val')
  y_mis_long$arm[y_mis_long$time=='t1'] <- 'trt'
  y_mis_long$val[y_mis_long$arm=='act'] <- NA
  y_mis_long <- na.omit(y_mis_long)
  
  m1 <- lme(val~time, 
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
  se_dj_mmrm <- est[2]
  low_dj_mmrm <- est[6]
  upp_dj_mmrm <- est[7]
  cover_dj_mmrm <- low_dj_mmrm<true_dj & true_dj<upp_dj_mmrm
  out <- c(out, theta_dj_mmrm, se_dj_mmrm, cover_dj_mmrm)
  
  
  # MI
  MI <- 50
  z <- matrix(c(1,1,0,1), nrow=2)
  sigma <- z%*%getVarCov(m1)%*%t(z) + diag(2)*m1$sigma^2
  sigma_mu <- z%*%m1$varFix%*%t(z)
  z_mu <- matrix(c(1,0,1,1),nrow=2,byrow=T)
  mu_act <- z_mu %*% summary(m1)$coef$fixed
  theta_dj_mi_j <- vartotal_dj_mi_j <- rep(NA, MI)
  for(j in 1:MI){
    mu_act_tilde <- rmvnorm(1, mu_act, sigma_mu)
    sigma_tilde <- riwish(N-1, sigma)*(N-1)

    y_mis$t2[ind_df] <- rnorm(length(ind_df),
                      mu_act_tilde[2]+sigma_tilde[1,2]/sigma_tilde[1,1]*(y_mis$t1[ind_df]-mu_act_tilde[1]),
                      sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))

    theta_dj_mi_j[j] <- mean(y_mis$t2)
    vartotal_dj_mi_j[j] <- var(y_mis$t2)/N
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

  
  
  ## De Facto missing
  # MMRM
  y_mis <- y_df
  ind_mis <- which(ind==3)
  y_mis$arm[ind_mis] <- 'mis'
  y_mis$t2[ind_mis] <- NA
  y_mis_long <- melt(y_mis,id.vars = c('id','arm'))
  colnames(y_mis_long) <- c('id','arm','time','val')
  y_mis_long$arm[y_mis_long$time=='t1'] <- 'act'
  y_mis_long$val[y_mis_long$arm=='trt'] <- NA
  y_mis_long <- na.omit(y_mis_long)
  
  m2 <- lme(val~time, 
            random=~1+time|id, 
            method='REML', 
            data=y_mis_long,
            control=lmeControl(returnObject=TRUE))
  est2 <- estimable(m2,
                   cm=cbind(
                     '(Intercept)'=1,
                     'timet2'=1
                   ),
                   conf.int=0.95)
  est2 <- as.numeric(est2)
  theta_df_mmrm <- est[1]*(1-r/N)+est2[1]*(r/N)
  se_df_mmrm <- sqrt(((N-r-1)*est[2]^2 + (r-1)*est2[2]^2)/(N-2))
  low_df_mmrm <- theta_df_mmrm - qt(1-alpha/2, N-2)*se_df_mmrm
  upp_df_mmrm <- theta_df_mmrm + qt(1-alpha/2, N-2)*se_df_mmrm
  cover_df_mmrm <- low_df_mmrm<true_df & true_df<upp_df_mmrm
  out <- c(out, theta_df_mmrm, se_df_mmrm, cover_df_mmrm)
  
  # MI
  MI <- 50
  z <- matrix(c(1,1,0,1), nrow=2)
  sigma <- z%*%getVarCov(m2)%*%t(z) + diag(2)*m2$sigma^2
  sigma_mu <- z%*%m2$varFix%*%t(z)
  z_mu <- matrix(c(1,0,1,1),nrow=2,byrow=T)
  mu_act <- z_mu %*% summary(m2)$coef$fixed
  theta_df_mi_j <- vartotal_df_mi_j <- rep(NA, MI)
  for(j in 1:MI){
    mu_act_tilde <- rmvnorm(1, mu_act, sigma_mu)
    sigma_tilde <- riwish(N-1, sigma)*(N-1)

    y_mis$t2[ind_mis] <- rnorm(length(ind_mis),
                      mu_act_tilde[2]+sigma_tilde[1,2]/sigma_tilde[1,1]*(y_mis$t1[ind_mis]-mu_act_tilde[1]),
                      sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))

    theta_df_mi_j[j] <- mean(y_mis$t2)
    vartotal_df_mi_j[j] <- var(y_mis$t2)/N
  }
  theta_df_mi <- mean(theta_df_mi_j)
  W <- mean(vartotal_df_mi_j)
  B <- var(theta_df_mi_j)
  se_df_mi <-  sqrt(W + (1+1/MI)*B)
  df <- (1+MI/(MI+1)*W/B)*(MI-1)
  low_df_mi <- theta_df_mi - qt(1-alpha/2, df)*se_df_mi
  upp_df_mi <- theta_df_mi + qt(1-alpha/2, df)*se_df_mi
  cover_df_mi <- low_df_mi<true_df & true_df<upp_df_mi

  out <- c(out, theta_df_mi, se_df_mi, cover_df_mi)
})

stopCluster(cl)

round(rowMeans(result),3)
round(sqrt(apply(result,1,var)),3)
