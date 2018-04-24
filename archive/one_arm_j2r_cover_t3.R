library(MASS)
library(nlme)
library(reshape2)
library(parallel)
library(invgamma)
library(MCMCpack)
library(Hmisc)
library(gmodels)

rm(list = ls())


mu_trt <- c(1, 2, 4)
mu_act1 <- c(.8, 1, 1)
mu_act2 <- c(.5, 2, 1)
mu_mis1 <- c(.3, 1, 1)
mu_mis2 <- c(.1, 2, 1)

var <- diag(c(2, 2.5, 3))
cor <- matrix(c(1, .7, .3,
                .7, 1, .5,
                .3, .5, 1), nrow=3)
vv <- var %*% cor %*% t(var)
N <- 100
alpha <- .05
# pi <- .7
# pi_c <- .15
# omega <- .7
# omega_c <- .15

pi <- .2
pi_a1 <- .2
pi_a2 <- .2
pi_m1 <- .2
pi_m2 <- .2


S <- 5000
true <- mu_trt[3]*pi + mu_act1[3]*pi_a1 + mu_act2[3]*pi_a2 + mu_mis1[3]*pi_m1 + mu_mis2[3]*pi_m2



no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl,varlist=c('S', 'N', 'vv', 'alpha', 'true',
                           'pi', 'pi_a1', 'pi_a2', 'pi_m1', 'pi_m2',
                           'mu_trt', 'mu_act1', 'mu_act2', 'mu_mis1', 'mu_mis2'))
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
  n1 <- sample(1:5,N,replace=TRUE,prob=c(pi,pi_a1,pi_a2,pi_m1,pi_m2))  # partition size
  n1_trt <- sum(n1==1)
  n1_act1 <- sum(n1==2)
  n1_act2 <- sum(n1==3)
  n1_mis1 <- sum(n1==4)
  n1_mis2 <- sum(n1==5)
  
  if(n1_trt*n1_act1*n1_act2*n1_mis1*n1_mis2 == 0){
    n1 <- sample(1:5,N,replace=TRUE,prob=c(pi,pi_a1,pi_a2,pi_m1,pi_m2))  # partition size
    n1_trt <- sum(n1==1)
    n1_act1 <- sum(n1==2)
    n1_act2 <- sum(n1==3)
    n1_mis1 <- sum(n1==4)
    n1_mis2 <- sum(n1==5)
  }
  
  
  y_trt <- y_act1 <- y_act1 <- y_mis1 <- y_mis2 <- data.frame()
  
  ## treatment on treatment ##
  y_trt <- data.frame(id=as.character(1:n1_trt), 
                      mvrnorm(n1_trt, mu_trt, vv), 
                      dose = 'treatment', group = 'one')
  colnames(y_trt) <- c('id', 't1', 't2', 't3', 'dose', 'group')
  
  
  ## treatment on active control at t1  and t2 ##
  y_act1 <- data.frame(id=as.character((n1_trt+1):(n1_trt+n1_act1)), 
                      mvrnorm(n1_act1, mu_act1, vv), 
                      dose = 'active1', group='one')
  colnames(y_act1) <- c('id', 't1', 't2', 't3', 'dose', 'group')
  
  
  ## treatment on active control at t2 ##
  y_act2 <- data.frame(id=as.character((n1_trt+n1_act1+1):(n1_trt+n1_act1+n1_act2)), 
                       mvrnorm(n1_act2, mu_act2, vv), 
                       dose = 'active2', group='one')
  colnames(y_act2) <- c('id', 't1', 't2', 't3', 'dose', 'group') 
  
  ## treatment missing at t2 and t2 ##
  y_mis1 <- data.frame(id=as.character((n1_trt+n1_act1+n1_act2+1):(n1_trt+n1_act1+n1_act2+n1_mis1)),
                      mvrnorm(n1_mis1, mu_mis1, vv), 
                      dose='miss1', group='one')
  colnames(y_mis1) <- c('id', 't1', 't2', 't3', 'dose', 'group') 
  
  ## treatment missing at t2 ##
  y_mis2 <- data.frame(id=as.character((n1_trt+n1_act1+n1_act2+n1_mis1+1):N),
                       mvrnorm(n1_mis2, mu_mis2, vv), 
                       dose='miss2', group='one')
  colnames(y_mis2) <- c('id', 't1', 't2', 't3', 'dose', 'group') 
  
  y_hyp <- rbind(y_trt, 
                 y_act1, y_act2,
                 y_mis1, y_mis2)  # The hypothetical complete data
  y_obs <- rbind(y_trt, 
                 y_act1, y_act2)   # The observed complete data
  y_obs2 <- rbind(y_trt, 
                  y_act1, y_act2,
                  y_mis2)   # The observed complete data at t2
  
  r <- n1_trt
  q1 <- n1_trt + n1_act1
  q2 <- n1_trt + n1_act1 + n1_act2
  
  # If completely observed
  theta <- mean(y_hyp$t3)
  vartotal <- var(y_hyp$t3)
  se <- sqrt(vartotal/N)
  low <- theta - qt(1-alpha/2, N-1)*se
  upp <- theta + qt(1-alpha/2, N-1)*se
  cover <- low<true & true<upp
  
  # If partial observed
  theta_obs <- mean(y_obs$t3)
  vartotal_obs <- var(y_obs$t3)
  se_obs <- sqrt(vartotal_obs/q2)
  low_obs <- theta_obs - qt(1-alpha/2, q2-1)*se_obs
  upp_obs <- theta_obs + qt(1-alpha/2, q2-1)*se_obs
  cover_obs <- low_obs<true & true<upp_obs
  
  ## MMRM ##
  y_obs_long <- melt(y_obs)
  colnames(y_obs_long) <- c('id','dose','group','time','val')
  
  m1 <- lme(val~time*dose, 
            random=~1+time|id, 
            method='REML', 
            data=y_obs_long,
            control=lmeControl(returnObject=TRUE))
  
  ## MI Rubin ##
  MI <- 20
  mu_act1 <- colMeans(y_act1[,2:4])
  mu_act2 <- colMeans(y_act2[,2:4])
  
  z <- matrix(c(1, 0, 0,
                1, 1, 0,
                1, 0, 1), nrow=3, byrow=T)
  sigma <- z%*%getVarCov(m1)%*%t(z) + diag(3)*m1$sigma^2
  
  theta_mi_j <- vartotal_mi_j <- varcond_mi_j <- rep(NA, MI)
  for(j in 1:MI){
    mu_act1_tilde <- mvrnorm(1, mu_act1, sigma/(n1_act1-1))
    mu_act2_tilde <- mvrnorm(1, mu_act2, sigma/(n1_act2-1))
    
    sigma_tilde <- riwish(q2-3, sigma)*(q2-3)
    
    # Impute missing at both t2 and t3
    sigma_tilde_12 <- t(as.matrix(sigma_tilde[1, 2:3]))
    sigma_tilde_21 <- as.matrix(sigma_tilde[2:3, 1])
    sigma_tilde_11 <- as.matrix(sigma_tilde[1, 1])
    sigma_tilde_22 <- as.matrix(sigma_tilde[2:3, 2:3])
    
    y_mis1[,3:4] <- t(apply(t(mu_act1_tilde[2:3]+
                            (sigma_tilde_21%*%solve(sigma_tilde_11)%*%(y_mis1$t1-mu_act1_tilde[1]))), 
                          1, function(x){
                            mvrnorm(1, x, sigma_tilde_22-sigma_tilde_21%*%solve(sigma_tilde_11)%*%sigma_tilde_12)}
    ))
    
    
    # Impute missing at t3
    sigma_tilde_12 <- as.matrix(sigma_tilde[1:2, 3])
    sigma_tilde_21 <- t(as.matrix(sigma_tilde[3, 1:2]))
    sigma_tilde_11 <- as.matrix(sigma_tilde[1:2, 1:2])
    sigma_tilde_22 <- as.matrix(sigma_tilde[3, 3])

    y_mis2[,4] <- rnorm(n1_mis2, 
                        mu_act2_tilde[3]+sigma_tilde_21%*%solve(sigma_tilde_11)%*%t(y_mis2[,2:3]-mu_act1_tilde[1:2]), 
                        sqrt(sigma_tilde_22-sigma_tilde_21%*%solve(sigma_tilde_11)%*%sigma_tilde_12))
    
    y3_mi <- c(y_obs$t3, y_mis1$t3, y_mis2$t3)
    theta_mi_j[j] <- mean(y3_mi)
    vartotal_mi_j[j] <- var(y3_mi)/N
  }
  theta_mi <- mean(theta_mi_j)
  W <- mean(vartotal_mi_j)
  B <- var(theta_mi_j)
  vartotal_mi <-  W + (1+1/MI)*B
  se_mi <- sqrt(vartotal_mi)
  df <- (1+MI/(MI+1)*W/B)*(MI-1)
  low_mi <- theta_mi - qt(1-alpha/2, df)*se_mi
  upp_mi <- theta_mi + qt(1-alpha/2, df)*se_mi
  cover_mi <- low_mi<true & true<upp_mi
  
  r <- c(theta, vartotal, se, cover,
         theta_obs, vartotal_obs, se_obs, cover_obs,
         theta_mi, vartotal_mi, se_mi, cover_mi)
})

stopCluster(cl)

rowMeans(result)
sqrt(apply(result,1,var) )
