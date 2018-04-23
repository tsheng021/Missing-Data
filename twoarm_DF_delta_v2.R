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
# psi2 <- c(-.45, -.4, 0)
# psi1p <- c(.4, -.4, 0)
# psi2p <- c(0, -.4, 0)

# S2
# psi1 <- c(.4, -.4, 0)
# psi2 <- c(-.45, -.4, 0)
# psi1p <- c(.4, -.4, 0)
# psi2p  <- c(.4, 0, -.4)

# S3
# psi1 <- c(1.6, 0, -.4)
# psi2 <- c(-.45, -.4, 0)
# psi1p <- c(1.6, 0, -.4)
# psi2p  <- c(0, -.4, 0)
# S2
psi1 <- c(1.6, 0, -.4)
psi2 <- c(-.45, -.4, 0)
psi1p <- c(1.6, 0, -.4)
psi2p  <- c(.4, 0, -.4)



px1 <- psi1 %*% c(1,mu_dj[2:3])
px2 <- psi2 %*% c(1,mu_df[2:3])
px1p <- psi1p %*% c(1,mu_pj[2:3])
px2p <- psi2p %*% c(1,mu_pf[2:3])
p1 <- 1/(1+exp(-px1))
p2 <- 1/(1+exp(-px2))
p1p <- 1/(1+exp(-px1p))
p2p <- 1/(1+exp(-px2p))
p1
p2
p1p
p2p
S <- 5000

# df <- c(NA,S)
# for(i in 1:S){
#   # Treatment arm
#   y_dj <- mvrnorm(N, mu_dj, vv)
#   y_dj <- data.frame(id=1:N, t0=y_dj[,1], t1=y_dj[,2], t2=y_dj[,3], arm='trt',stringsAsFactors=FALSE)
#   px <- data.frame(p1=as.matrix(cbind(1,y_dj[,3:4])) %*% t(t(psi1)),
#                    p2=as.matrix(cbind(1,y_dj[,3:4])) %*% t(t(psi2)))
#   pt <- 1/(1+exp(-px))
#   
#   ind <- rbinom(N, 1, pt[,1])+1
#   ind[ind==2] <- ind[ind==2] + rbinom(sum(ind==2), 1, pt[,2])
#   y_df <- y_dj
#   
#   ind_df <- which(ind==2 | ind==3)
#   y_df$t2[ind_df] <- rnorm(length(ind_df), 
#                            mean=mu_df[3]+vv[3,1:2]%*%solve(vv[1:2,1:2])%*%t(y_df[ind_df,2:3]-mu_dj[1:2]),
#                            sd=sqrt(vv[3,3]-vv[3,1:2]%*%solve(vv[1:2,1:2])%*%vv[1:2,3]))
#   
#   # Placebo arm
#   y_pj <- mvrnorm(N, mu_pj, vv)
#   y_pj <- data.frame(id=1:N, t0=y_pj[,1], t1=y_pj[,2], t2=y_pj[,3], arm='pla',stringsAsFactors=FALSE)
#   px <- data.frame(p1=as.matrix(cbind(1,y_pj[,3:4])) %*% t(t(psi1p)),
#                    p2=as.matrix(cbind(1,y_pj[,3:4])) %*% t(t(psi2p)))
#   pt <- 1/(1+exp(-px))
#   ind <- rbinom(N, 1, pt[,1])+1
#   ind[ind==2] <- ind[ind==2] + rbinom(sum(ind==2), 1, pt[,2])
#   y_pf <- y_pj
#   
#   ind_df <- which(ind==2 | ind==3)
#   y_pf$t2[ind_df] <- rnorm(length(ind_df), 
#                            mean=mu_pf[3]+vv[3,1:2]%*%solve(vv[1:2,1:2])%*%t(y_pf[ind_df,2:3]-mu_pj[1:2]),
#                            sd=sqrt(vv[3,3]-vv[3,1:2]%*%solve(vv[1:2,1:2])%*%vv[1:2,3]))
#   
#   
#   df[i] <- mean(y_df$t2)-mean(y_pf$t2)
#   
# }
# 
true_dj <- mu_dj[2]-mu_pj[2]
# true_df <- mean(df)
true_df <- 0




no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl,varlist=c('S', 'N', 'vv', 'alpha', 
                           'true_dj','true_df',
                           'psi1', 'psi2', 
                           'psi1p', 'psi2p',
                           'mu_pj','mu_pf',
                           'mu_dj','mu_df'))
clusterSetRNGStream(cl, 1234)
clusterEvalQ(cl, {
  library(MASS)
  library(nlme)
  library(reshape2)
})


result <- parSapply(cl, 1:S, function(x){
  out <- c()
  y_dj <- mvrnorm(N, mu_dj, vv)
  y_dj <- data.frame(id=1:N, t0=y_dj[,1], t1=y_dj[,2], t2=y_dj[,3], arm='trt',stringsAsFactors=FALSE)
  px <- as.matrix(cbind(1,y_dj[,3:4])) %*% t(t(psi1))
  pt <- 1/(1+exp(-px))
  ind <- rbinom(N, 1, pt)+1
  if (sum(ind==2)==0) ind <- rbinom(N, 1, pt)+1
  ind_df <- which(ind==2)
  y_df <- y_dj
  y_df$arm <- 'trt'
  y_df$arm[ind_df] <- 'act'
  y_df$t2[ind_df] <- rnorm(length(ind_df), 
                           mean=mu_df[3]+vv[3,1:2]%*%solve(vv[1:2,1:2])%*%t(sweep(y_df[ind_df,2:3], 2, mu_dj[1:2])),
                           sd=sqrt(vv[3,3]-vv[3,1:2]%*%solve(vv[1:2,1:2])%*%vv[1:2,3]))
  px <- as.matrix(cbind(1,y_dj[ind_df,3:4])) %*% t(t(psi2))
  pt <- 1/(1+exp(-px))
  ind[ind==2] <- ind[ind==2] + rbinom(sum(ind==2), 1, pt)
  if (sum(ind==3)==0) ind[ind==2]<-ind[ind==2] + rbinom(sum(ind==2), 1, pt)
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
  ind <- rbinom(N, 1, pt)+1
  if (sum(ind==2)==0) ind <- rbinom(N, 1, pt)+1
  y_pf <- y_pj
  indp_df <- which(ind==2)
  y_pf$arm <- 'pla'
  y_pf$arm[indp_df] <- 'act'
  y_pf$t2[indp_df] <- rnorm(length(indp_df), 
                            mean=mu_pf[3]+vv[3,1:2]%*%solve(vv[1:2,1:2])%*%t(sweep(y_pf[indp_df,2:3], 2, mu_pj[1:2])),
                            sd=sqrt(vv[3,3]-vv[3,1:2]%*%solve(vv[1:2,1:2])%*%vv[1:2,3]))
  px <- as.matrix(cbind(1,y_pj[indp_df,3:4])) %*% t(t(psi2p))
  pt <- 1/(1+exp(-px))
  ind[ind==2] <- ind[ind==2] + rbinom(sum(ind==2), 1, pt)
  if (sum(ind==3)==0) ind[ind==2]<-ind[ind==2] + rbinom(sum(ind==2), 1, pt)
  y_misp <- y_pf
  indp_mis <- which(ind==3)
  y_misp$t2[indp_mis] <- NA
  y_misp_long <- melt(y_misp,id.vars = c('id','arm','t0'))
  colnames(y_misp_long) <- c('id','arm','bl','time','val')
  y_misp_long <- na.omit(y_misp_long)
  
  l <- N-length(indp_df)
  m <- N-l-length(indp_mis)
  
  # Total dataset
  tot_mis_long <- rbind(data.frame(y_mis_long,group='treatment'), data.frame(y_misp_long,group='placebo')) 
  tot_mis_long$d1 <- as.numeric(tot_mis_long$arm=='act' & 
                                  tot_mis_long$group=='treatment' &
                                  tot_mis_long$time=='t2')
  tot_mis_long$d2 <- as.numeric(tot_mis_long$arm=='act' & 
                                  tot_mis_long$group=='placebo' &
                                  tot_mis_long$time=='t2')
  
  # MMRM
  m1 <- lme(val~bl+time+group+group*time+bl*time+d1+d2, 
            random=~(-1)+time|id, 
            method='REML', 
            data=tot_mis_long,
            control=lmeControl(returnObject=TRUE))
  bl_bar <- mean(c(y_mis$t0,y_misp$t0),na.rm=T)
  
  # MAR-MMRM
  coeff <- summary(m1)$coef$fixed
  D_mar <- c(0, 0, 0, -1, 1-r/N, -1+l/N, -1, 0)
  theta_mar_mmrm <- D_mar%*%coeff
  se_mmrm <- sqrt(D_mar%*%vcov(m1)%*%D_mar)
  se_mar_mmrm <- sqrt(se_mmrm^2 + 
                        coeff[5]^2*r*(N-r)/N^3+
                        coeff[6]^2*l*(N-l)/N^3)
  low_mar_mmrm <- theta_mar_mmrm - qt(1-alpha/2, N-2)*se_mar_mmrm
  upp_mar_mmrm <- theta_mar_mmrm + qt(1-alpha/2, N-2)*se_mar_mmrm
  cover_mar_mmrm <- low_mar_mmrm<true_df & true_df<upp_mar_mmrm
  t <- abs(theta_mar_mmrm/se_mar_mmrm)
  pval_mar_mmrm <- ((1-pt(t,2*N))<.025)*100
  out <- c(out, theta_mar_mmrm, se_mar_mmrm, cover_mar_mmrm, pval_mar_mmrm)
  
  # MI
  MI <- 20
  sigma <- getVarCov(m1) + diag(2)*m1$sigma^2
  theta_mar_mi_j <- vartotal_mar_mi_j <- rep(NA, MI)
  for(j in 1:MI){
    y_mis$t2[ind_mis] <- NA
    y_misp$t2[indp_mis] <- NA
    for(k in 1:N){
      if(is.na(y_mis$t2[k])){
        z_act <- matrix(c(1,y_mis$t0[k],0,0,0,0,0,y_mis$t0[k],
                          1,y_mis$t0[k],1,0,1,0,0,y_mis$t0[k]),nrow=2,byrow=T)
        mu_act <- z_act %*% coeff
        mu_act_tilde <- mvrnorm(1, mu_act, sigma/(N-r-1))
        sigma_tilde <- matrix(rWishart(1, 2*N-l-r-2, sigma),nrow=2)/(2*N-l-r-2)
        y_mis$t2[k] <- rnorm(1,
                             mu_act_tilde[2]+sigma_tilde[2,1]/sigma_tilde[1,1]*(y_mis$t1[k]-mu_act_tilde[1]),
                             sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
        
      }
      if(is.na(y_misp$t2[k])){
        z_act <- matrix(c(1,y_misp$t0[k],0,1,0,0,0,y_misp$t0[k],
                          1,y_misp$t0[k],1,1,0,1,1,y_misp$t0[k]),nrow=2,byrow=T)
        
        mu_act <- z_act %*% coeff
        mu_act_tilde <- mvrnorm(1, mu_act, sigma/(N-l-1))
        sigma_tilde <- matrix(rWishart(1, 2*N-l-r-2, sigma),nrow=2)/(2*N-l-r-2)
        y_misp$t2[k] <- rnorm(1,
                             mu_act_tilde[2]+sigma_tilde[2,1]/sigma_tilde[1,1]*(y_misp$t1[k]-mu_act_tilde[1]),
                             sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
        
      }
    }
    theta_mar_mi_j[j] <- mean(y_mis$t2)-mean(y_misp$t2)
    vartotal_mar_mi_j[j] <- (var(y_mis$t2)+var(y_misp$t2))/N
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
  t <- abs(theta_mar_mi/se_mar_mi)
  pval_mar_mi <- (1-pt(t, df)<.025)*100
  out <- c(out, theta_mar_mi, se_mar_mi, cover_mar_mi, pval_mar_mi)
})

stopCluster(cl)

round(rowMeans(result),2)
round(sqrt(apply(result,1,var)),2)
