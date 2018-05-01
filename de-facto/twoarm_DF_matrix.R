library(MASS)
library(nlme)
library(reshape2)
library(parallel)
library(emmeans)

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
psi1 <- c(.4, -.4, 0)
psi2 <- c(-.3, -.4, 0)
psi1p <- c(-.4, .4, 0)
psi2p <- c(-1.6, .4, 0)

S <- 1000

df <- pval <-
  switch_trt <- miss_trt <-
  switch_pla <- miss_pla <-
  c(NA,S)
for(i in 1:S){
  # Treatment arm
  y_dj <- mvrnorm(N, mu_dj, vv)
  y_dj <- data.frame(id=1:N, t0=y_dj[,1], t1=y_dj[,2], t2=y_dj[,3],stringsAsFactors=FALSE)
  px <- as.matrix(cbind(1,y_dj[,3:4])) %*% t(t(psi1))
  pt <- 1/(1+exp(-px))
  ind <- rbinom(N, 1, pt)+1
  while (sum(ind==2)==0) ind <- rbinom(N, 1, pt)+1
  ind_df <- which(ind==2)
  y_df <- y_dj
  y_df$g <- 1
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
  y_pj <- data.frame(id=1:N+N, t0=y_pj[,1], t1=y_pj[,2], t2=y_pj[,3],stringsAsFactors=FALSE)
  px <- as.matrix(cbind(1,y_pj[,3:4])) %*% t(t(psi1p))
  pt <- 1/(1+exp(-px))
  indp <- rbinom(N, 1, pt)+1
  while (sum(indp==2)==0) indp <- rbinom(N, 1, pt)+1
  indp_df <- which(indp==2)
  y_pf <- y_pj
  y_pf$g <- 0
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
  
  tot_f <- rbind(y_df, y_pf)
  tot_f$s <- c(ind, indp)
  tot_long <- melt(tot_f,id.var=c('id','t0','g','s'))
  tot_long$g <- as.factor(tot_long$g)
  colnames(tot_long) <- c('id', 'bl', 'g', 's', 'time', 'val')
  # Create coding variable
  tot_long$s1g1 <- 0
  tot_long$s1g1[which(tot_long$s==1 & tot_long$g==1)] <- 1
  
  tot_long$s2g0 <- 0
  tot_long$s2g0[which(tot_long$s>1 & tot_long$g==0)] <- 1
  
  tot_long$s2g1 <- 0
  tot_long$s2g1[which(tot_long$s>1 & tot_long$g==1)] <- 1
  
  tot_long$s23g1 <- 0
  tot_long$s23g1[which(tot_long$s>1 & tot_long$g==1)] <- 1
  
  tot_long$s23g0 <- 0
  tot_long$s23g0[which(tot_long$s>1 & tot_long$g==0)] <- 1
  
  tot_long$s1g1t2 <- 0
  tot_long$s1g1t2[which(tot_long$s==1 & 
                          tot_long$g==1 & 
                          tot_long$time=='t2')] <- 1
  
  tot_long$s2g0t2 <- 0
  tot_long$s2g0t2[which(tot_long$s==2 & 
                          tot_long$g==0 & 
                          tot_long$time=='t2')] <- 1
  
  tot_long$s2g1t2 <- 0
  tot_long$s2g1t2[which(tot_long$s==2 & 
                          tot_long$g==1 & 
                          tot_long$time=='t2')] <- 1
  
  tot_long$s3g0t2 <- 0
  tot_long$s3g0t2[which(tot_long$s==3 & 
                          tot_long$g==0 & 
                          tot_long$time=='t2')] <- 1
  
  tot_long$s3g1t2 <- 0
  tot_long$s3g1t2[which(tot_long$s==3 & 
                          tot_long$g==1 & 
                          tot_long$time=='t2')] <- 1
  
  # MMRM
  m1 <- lme(val~bl+time+s1g1+s23g0+s23g1+bl:time+s1g1t2+s2g0t2+s2g1t2+s3g0t2+s3g1t2,
            random=~(-1)+time|id,
            method='REML',
            data=tot_long,
            control=lmeControl(returnObject=TRUE))
  
  m1_coef <- summary(m1)$coef$fixed[c(1,3,4,5,6,7,8,9,10,11)]
  sig_b <- vcov(m1)[c(1,3),c(1,3)]
  sig_g <- vcov(m1)[c(4,5,6,7,8,9,10,11),c(4,5,6,7,8,9,10,11)]
  cov_bg <- vcov(m1)[c(4,5,6,7,8,9,10,11), c(1,3)]
  
  l <- length(which(indp==1))
  m <- length(which(indp==2))
  n <- length(which(indp==3))
  r <- length(which(ind==1))
  q <- length(which(indp==2))
  p <- length(which(indp==3))
  
  sig_o <- matrix(c(m*(N-m)/N^3, -n*m/N^3, 0, 0, 
                    -n*m/N^3, n*(N-n)/N^3, 0, 0,
                    0, 0, q*(N-q)/N^3, -q*p/N^3,
                    0, 0, -q*p/N^3, p*(N-p)/N^3), nrow=4, byrow=T)
  
  ws <- length(which(indp>1))/N
  wm <- length(indp_mis)/length(indp_df)
  ps <- length(which(ind>1))/N
  pm <- length(ind_mis)/length(ind_df)
  
  U_marg <- matrix(c(1, 0, 0, 0,
                     1, 0, 1, 0,
                     1, 1, 0, 0,
                     1, 1, 1, 1), nrow=4, byrow=T)
  U_c1 <- matrix(c(1, 0,
                   1, 1, 
                   1, 0, 
                   1, 1), nrow=4, byrow=T)
  U_c2 <-  matrix(c(0, ws, 0, 0, 0, 0, 0, 0,
                    0, ws, 0, 0, (1-wm)*ws, 0, wm*ws, 0,
                    1-ps, 0, ps, 0, 0, 0, 0, 0,
                    1-ps, 0, ps, 1-ps, 0, (1-pm)*ps, 0, pm*ps), nrow=4, byrow=T)
  beta <- matrix(m1_coef[c(1,2)],ncol=1)
  gamma <- matrix(m1_coef[3:10], ncol=1)
  a1 <- solve(t(U_marg)%*%U_marg)%*%t(U_marg)%*%(U_c1%*%beta+U_c2%*%gamma)
  
  v_a1 <- solve(t(U_marg)%*%U_marg)%*%t(U_marg)
  o <- matrix(c(sum(indp==2), sum(indp==3), sum(ind==2), sum(ind==3)),
             ncol=1)
  X <- U_c2%*%gamma%*%t(o)%*%ginv(o%*%t(o))
  
  v_a1 <- c(0,1,0,1)%*%solve(t(U_marg)%*%U_marg)%*%t(U_marg)%*%
    (X%*%sig_o%*%t(X)+
       U_c1%*%sig_b%*%t(U_c1)+
       U_c2%*%sig_g%*%t(U_c2)+
       U_c1%*%t(cov_bg)%*%t(U_c2)+
       U_c2%*%cov_bg%*%t(U_c1))%*%
    t(solve(t(U_marg)%*%U_marg)%*%t(U_marg))%*%t(t(c(0,1,0,1)))
  
  df[i] <- a1[2] + a1[4]
  pval[i] <- ((1-pt(df[i]/sqrt(v_a1),2*N))<.025)*100
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


true_df <- 0

no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl,varlist=c('S', 'N', 'vv', 'alpha', 
                           'true_df',
                           'psi1', 'psi2', 
                           'psi1p', 'psi2p',
                           'mu_pj','mu_pf',
                           'mu_dj','mu_df'))
clusterSetRNGStream(cl, 4307)
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
  y_dj <- data.frame(id=1:N, t0=y_dj[,1], t1=y_dj[,2], t2=y_dj[,3], g='trt',stringsAsFactors=FALSE)
  px <- as.matrix(cbind(1,y_dj[,3:4])) %*% t(t(psi1))
  pt <- 1/(1+exp(-px))
  ind <- rbinom(N, 1, pt)+1
  while (sum(ind==2)==0) ind <- rbinom(N, 1, pt)+1
  ind_df <- which(ind==2)
  y_df <- y_dj
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
  ind_mis <- which(ind==3)
  y_df$s <- ind
  y_mis <- y_df
  y_mis$t2[ind_mis] <- NA
  y_mis_long <- melt(y_mis,id.vars = c('id','g','t0','s'))
  colnames(y_mis_long) <- c('id','g','bl','s','time','val')
  y_mis_long <- na.omit(y_mis_long)
  
  r <- N-length(ind_df)
  q <- N-r-length(ind_mis)
  
  
  # Placebo arm
  y_pj <- mvrnorm(N, mu_pj, vv)
  y_pj <- data.frame(id=1:N+N, t0=y_pj[,1], t1=y_pj[,2], t2=y_pj[,3], g='pla',stringsAsFactors=FALSE)
  px <- as.matrix(cbind(1,y_pj[,3:4])) %*% t(t(psi1p))
  pt <- 1/(1+exp(-px))
  indp <- rbinom(N, 1, pt)+1
  while (sum(indp==2)==0) indp <- rbinom(N, 1, pt)+1
  y_pf <- y_pj
  indp_df <- which(indp==2)
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
  indp_mis <- which(indp==3)
  y_pf$s <- indp
  y_misp <- y_pf
  y_misp$t2[indp_mis] <- NA
  y_misp_long <- melt(y_misp,id.vars = c('id','g','t0','s'))
  colnames(y_misp_long) <- c('id','g','bl','s','time','val')
  y_misp_long <- na.omit(y_misp_long)
  
  l <- N-length(indp_df)
  m <- N-l-length(indp_mis)
  
  # Total dataset
  tot_mis_long <- rbind(y_mis_long, y_misp_long)
  
  # Create coding variable
  tot_mis_long$s1g1 <- 0
  tot_mis_long$s1g1[which(tot_mis_long$s==1 & tot_mis_long$g=='trt')] <- 1
  
  tot_mis_long$s2g0 <- 0
  tot_mis_long$s2g0[which(tot_mis_long$s==2 & tot_mis_long$g=='pla')] <- 1
  
  tot_mis_long$s2g1 <- 0
  tot_mis_long$s2g1[which(tot_mis_long$s==2 & tot_mis_long$g=='trt')] <- 1
  
  tot_mis_long$s23g1 <- 0
  tot_mis_long$s23g1[which(tot_mis_long$s>1 & tot_mis_long$g=='trt')] <- 1
  
  tot_mis_long$s23g0 <- 0
  tot_mis_long$s23g0[which(tot_mis_long$s>1 & tot_mis_long$g=='pla')] <- 1
  
  tot_mis_long$s1g1t2 <- 0
  tot_mis_long$s1g1t2[which(tot_mis_long$s==1 & 
                              tot_mis_long$g=='trt' & 
                              tot_mis_long$time=='t2')] <- 1
  
  tot_mis_long$s2g0t2 <- 0
  tot_mis_long$s2g0t2[which(tot_mis_long$s==2 & 
                              tot_mis_long$g=='pla' & 
                              tot_mis_long$time=='t2')] <- 1
  
  tot_mis_long$s2g1t2 <- 0
  tot_mis_long$s2g1t2[which(tot_mis_long$s==2 & 
                              tot_mis_long$g=='trt' & 
                              tot_mis_long$time=='t2')] <- 1
  
  # MMRM
  m1 <- lme(val~bl+time+s1g1+s23g0+s23g1+bl:time+s1g1t2+s2g0t2+s2g1t2,
            random=~(-1)+time|id,
            method='REML',
            data=tot_mis_long,
            control=lmeControl(returnObject=TRUE))
  # bl_bar <- mean(c(y_df$t0,y_pf$t0))
  bl_bar <- mean(y_df$t0[ind_mis])
  
  # MAR-MMRM
  b <- summary(m1)$coef$fixed
  D_mar <- c(0, 0, 0, r/N, -(N-l)/N, (N-r)/N, r/N, -(N-l)/N, (N-r)/N, 0)
  theta_mar_mmrm <- D_mar%*%b
  se_mmrm <- sqrt(D_mar%*%vcov(m1)%*%D_mar)
  se_mar_mmrm <- sqrt(se_mmrm^2 +
                        r*(N-r)/N^3*(b[4]+b[7]-b[6]-b[9])^2+
                        l*(N-l)/N^3*(b[5]+b[8])^2)
  low_mar_mmrm <- theta_mar_mmrm - qt(1-alpha/2, N-2)*se_mar_mmrm
  upp_mar_mmrm <- theta_mar_mmrm + qt(1-alpha/2, N-2)*se_mar_mmrm
  cover_mar_mmrm <- low_mar_mmrm<true_df & true_df<upp_mar_mmrm
  t <- theta_mar_mmrm/se_mar_mmrm
  pval_mar_mmrm <- ((1-pt(t,2*N))<.025)*100
  out <- c(out, theta_mar_mmrm, se_mar_mmrm, cover_mar_mmrm, pval_mar_mmrm)
  
  # J2R-MMRM
  b <- summary(m1)$coef$fixed
  D_j2r <- c(0, 0, 0, r/N, (l-r-q)/N, q/N, (N-r)/N, (l-r-q)/N, q/N, 0)
  theta_j2r_mmrm <- D_j2r%*%b
  se_mmrm <- sqrt(D_j2r%*%vcov(m1)%*%D_j2r)
  se_j2r_mmrm <- sqrt(se_mmrm^2 +
                        r*(N-r)/N^3*(b[4]+b[7]-b[5]-b[8])^2+
                        q*(N-q)/N^3*(b[6]+b[9]-b[5]-b[8])^2-
                        2*r*q/N^3*(b[4]+b[7]-b[5]-b[8])*(b[6]+b[9]-b[5]-b[8])+
                        l*(N-l)/N^3*(b[5]+b[8])^2)
  low_j2r_mmrm <- theta_j2r_mmrm - qt(1-alpha/2, N-2)*se_j2r_mmrm
  upp_j2r_mmrm <- theta_j2r_mmrm + qt(1-alpha/2, N-2)*se_j2r_mmrm
  cover_j2r_mmrm <- low_j2r_mmrm<true_df & true_df<upp_j2r_mmrm
  t <- theta_j2r_mmrm/se_j2r_mmrm
  pval_j2r_mmrm <- ((1-pt(t,2*N))<.025)*100
  out <- c(out, theta_j2r_mmrm, se_j2r_mmrm, cover_j2r_mmrm, pval_j2r_mmrm)
  
  # J2B-MMRM
  b <- summary(m1)$coef$fixed
  D_j2b <- c(-(N-r-q)/N, -(N-r-q)/N*bl_bar, -(N-r-q)/N, r/N, -(N-l)/N, q/N, r/N, -(N-l)/N, q/N,-(N-r-q)/N*bl_bar)
  theta_j2b_mmrm <- D_j2b%*%b+(N-r-q)/N*bl_bar
  se_mmrm <- sqrt(D_j2b%*%vcov(m1)%*%D_j2b)
  se_j2b_mmrm <- sqrt(se_mmrm^2 +
                        r*(N-r)/N^3*(b[1]+b[3]-b[4]-b[7]+(1+b[2]+b[10])*bl_bar)^2+
                        q*(N-q)/N^3*(b[1]+b[3]-b[6]-b[9]+(1+b[2]+b[10])*bl_bar)^2-
                        2*r*q/N^3*(b[1]+b[3]-b[4]-b[7]+(1+b[2]+b[10])*bl_bar)*(b[1]+b[3]-b[6]-b[9]+(1+b[2]+b[10])*bl_bar)+
                        l*(N-l)/N^3*(b[5]+b[8])^2)
  low_j2b_mmrm <- theta_j2b_mmrm - qt(1-alpha/2, N-2)*se_j2b_mmrm
  upp_j2b_mmrm <- theta_j2b_mmrm + qt(1-alpha/2, N-2)*se_j2b_mmrm
  cover_j2b_mmrm <- low_j2b_mmrm<true_df & true_df<upp_j2b_mmrm
  t <- theta_j2b_mmrm/se_j2b_mmrm
  pval_j2b_mmrm <- ((1-pt(t,2*N))<.025)*100
  out <- c(out, theta_j2b_mmrm, se_j2b_mmrm, cover_j2b_mmrm, pval_j2b_mmrm)
  
  
  # SWP-MMRM
  b <- summary(m1)$coef$fixed
  D_swp <- c(0,0,0,r/N,(N-r-q-m)/N,(q-N+l+m)/N,r/N,(N-r-q-m)/N,(q-N+l+m)/N,0)
  theta_swp_mmrm <- D_swp%*%b
  se_mmrm <- sqrt(D_swp%*%vcov(m1)%*%D_swp)
  se_swp_mmrm <- sqrt(se_mmrm^2 +
                        r*(N-r)/N^3*(b[4]+b[7]-b[5]-b[8])^2+
                        q*(N-q)/N^3*(b[6]+b[9]-b[5]-b[8])^2-
                        2*r*q/N^3*(b[4]+b[7]-b[5]-b[8])*(b[6]+b[9]-b[5]-b[8])+
                        l*(N-l)/N^3*(b[6]+b[9])^2+m*(N-m)/N^3*(b[6]+b[9]-b[5]-b[8])^2-
                        2*l*m/N^3*(b[6]+b[9])*(b[6]+b[9]-b[5]-b[8]))
  low_swp_mmrm <- theta_swp_mmrm - qt(1-alpha/2, N-2)*se_swp_mmrm
  upp_swp_mmrm <- theta_swp_mmrm + qt(1-alpha/2, N-2)*se_swp_mmrm
  cover_swp_mmrm <- low_swp_mmrm<true_df & true_df<upp_swp_mmrm
  t <- theta_swp_mmrm/se_swp_mmrm
  pval_swp_mmrm <- ((1-pt(t,2*N))<.025)*100
  out <- c(out, theta_swp_mmrm, se_swp_mmrm, cover_swp_mmrm, pval_swp_mmrm)
  
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
        z_act <- matrix(c(1,y_mis$t0[k],0,0,0,1,0,0,0,0,
                          1,y_mis$t0[k],1,0,0,1,0,0,1,y_mis$t0[k]),nrow=2,byrow=T)
        mu_act <- z_act %*% b
        mu_act_tilde <- mvrnorm(1, mu_act, sigma)
        y_mis$t2[k] <- rnorm(1,
                             mu_act_tilde[2]+sigma_tilde[2,1]/sigma_tilde[1,1]*(y_mis$t1[k]-mu_act_tilde[1]),
                             sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
      }
      if(is.na(y_misp$t2[k])){
        z_act <- matrix(c(1,y_misp$t0[k],0,0,1,0,0,0,0,0,
                          1,y_misp$t0[k],1,0,1,0,0,1,0,y_misp$t0[k]),nrow=2,byrow=T)
        mu_act <- z_act %*% b
        mu_act_tilde <- mvrnorm(1, mu_act, sigma)
        y_misp$t2[k] <- rnorm(1,
                              mu_act_tilde[2]+sigma_tilde[2,1]/sigma_tilde[1,1]*(y_misp$t1[k]-mu_act_tilde[1]),
                              sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
        
      }
    }
    y_mi_tot <- rbind(y_mis,y_misp)
    y_mi <- melt(y_mi_tot,id.vars = c('id','g','t0','s'))
    colnames(y_mi) <- c('id','g','bl','s','time','val')
    lmi <- lme(val~bl+time+g+g:time+bl:time,
               random=~(-1)+time|id,
               method='REML',
               data=y_mi,
               control=lmeControl(returnObject=TRUE))
    lsm <- summary(lsmeans(lmi, list(pairwise ~ g:time), adjust='none', data=y_mi)[[2]][6])
    theta_mar_mi_j[j] <- as.numeric(-lsm[2])
    vartotal_mar_mi_j[j] <- as.numeric(lsm[3]^2)
    
    # Mild
    y_mis$t2[ind_mis] <- NA
    y_misp$t2[indp_mis] <- NA
    for(k in 1:N){
      if(is.na(y_mis$t2[k])){
        z_act <- matrix(c(1,y_mis$t0[k],0,0,0,1,0,0,0,0,
                          1,y_mis$t0[k],1,0,1,0,0,1,0,y_mis$t0[k]),nrow=2,byrow=T)
        mu_act <- z_act %*% b
        mu_act_tilde <- mvrnorm(1, mu_act, sigma)
        y_mis$t2[k] <- rnorm(1,
                             mu_act_tilde[2]+sigma_tilde[2,1]/sigma_tilde[1,1]*(y_mis$t1[k]-mu_act_tilde[1]),
                             sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
      }
      if(is.na(y_misp$t2[k])){
        z_act <- matrix(c(1,y_misp$t0[k],0,0,1,0,0,0,0,0,
                          1,y_misp$t0[k],1,0,1,0,0,1,0,y_misp$t0[k]),nrow=2,byrow=T)
        mu_act <- z_act %*% b
        mu_act_tilde <- mvrnorm(1, mu_act, sigma)
        y_misp$t2[k] <- rnorm(1,
                              mu_act_tilde[2]+sigma_tilde[2,1]/sigma_tilde[1,1]*(y_misp$t1[k]-mu_act_tilde[1]),
                              sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
      }
    }
    y_mi_tot <- rbind(y_mis,y_misp)
    y_mi <- melt(y_mi_tot,id.vars = c('id','g','t0','s'))
    colnames(y_mi) <- c('id','g','bl','s','time','val')
    lmi <- lme(val~bl+time+g+g:time+bl:time,
               random=~(-1)+time|id,
               method='REML',
               data=y_mi,
               control=lmeControl(returnObject=TRUE))
    lsm <- summary(lsmeans(lmi, list(pairwise ~ g:time), adjust='none', data=y_mi)[[2]][6])
    theta_mild_mi_j[j] <- as.numeric(-lsm[2])
    vartotal_mild_mi_j[j] <- as.numeric(lsm[3]^2)
    
    # Moderate
    y_mis$t2[ind_mis] <- NA
    y_misp$t2[indp_mis] <- NA
    for(k in 1:N){
      if(is.na(y_mis$t2[k])){
        z_act <- matrix(c(1,y_mis$t0[k],0,0,0,1,0,0,0,0,
                          1,y_mis$t0[k],1,0,1,0,0,1,0,y_mis$t0[k]),nrow=2,byrow=T)
        mu_act <- z_act %*% b
        mu_act[2] <- y_mis$t0[k]
        mu_act_tilde <- mvrnorm(1, mu_act, sigma)
        y_mis$t2[k] <- rnorm(1,
                             mu_act_tilde[2]+sigma_tilde[2,1]/sigma_tilde[1,1]*(y_mis$t1[k]-mu_act_tilde[1]),
                             sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
      }
      if(is.na(y_misp$t2[k])){
        z_act <- matrix(c(1,y_misp$t0[k],0,0,1,0,0,0,0,0,
                          1,y_misp$t0[k],1,0,1,0,0,1,0,y_misp$t0[k]),nrow=2,byrow=T)
        mu_act <- z_act %*% b
        mu_act_tilde <- mvrnorm(1, mu_act, sigma)
        y_misp$t2[k] <- rnorm(1,
                              mu_act_tilde[2]+sigma_tilde[2,1]/sigma_tilde[1,1]*(y_misp$t1[k]-mu_act_tilde[1]),
                              sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
      }
    }
    y_mi_tot <- rbind(y_mis,y_misp)
    y_mi <- melt(y_mi_tot,id.vars = c('id','g','t0','s'))
    colnames(y_mi) <- c('id','g','bl','s','time','val')
    lmi <- lme(val~bl+time+g+g:time+bl:time,
               random=~(-1)+time|id,
               method='REML',
               data=y_mi,
               control=lmeControl(returnObject=TRUE))
    lsm <- summary(lsmeans(lmi, list(pairwise ~ g:time), adjust='none', data=y_mi)[[2]][6])
    theta_moderate_mi_j[j] <- as.numeric(-lsm[2])
    vartotal_moderate_mi_j[j] <- as.numeric(lsm[3]^2)
    
    # Extreme
    y_mis$t2[ind_mis] <- NA
    y_misp$t2[indp_mis] <- NA
    for(k in 1:N){
      if(is.na(y_mis$t2[k])){
        z_act <- matrix(c(1,y_mis$t0[k],0,0,0,1,0,0,0,0,
                          1,y_mis$t0[k],1,0,1,0,0,1,0,y_mis$t0[k]),nrow=2,byrow=T)
        mu_act <- z_act %*% b
        mu_act_tilde <- mvrnorm(1, mu_act, sigma)
        y_mis$t2[k] <- rnorm(1,
                             mu_act_tilde[2]+sigma_tilde[2,1]/sigma_tilde[1,1]*(y_mis$t1[k]-mu_act_tilde[1]),
                             sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
      }
      if(is.na(y_misp$t2[k])){
        z_act <- matrix(c(1,y_misp$t0[k],0,0,1,0,0,0,0,0,
                          1,y_misp$t0[k],1,0,0,1,0,0,1,y_misp$t0[k]),nrow=2,byrow=T)
        mu_act <- z_act %*% b
        mu_act_tilde <- mvrnorm(1, mu_act, sigma)
        y_misp$t2[k] <- rnorm(1,
                              mu_act_tilde[2]+sigma_tilde[2,1]/sigma_tilde[1,1]*(y_misp$t1[k]-mu_act_tilde[1]),
                              sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
        
      }
    }
    y_mi_tot <- rbind(y_mis,y_misp)
    y_mi <- melt(y_mi_tot,id.vars = c('id','g','t0','s'))
    colnames(y_mi) <- c('id','g','bl','s','time','val')
    lmi <- lme(val~bl+time+g+g:time+bl:time,
               random=~(-1)+time|id,
               method='REML',
               data=y_mi,
               control=lmeControl(returnObject=TRUE))
    lsm <- summary(lsmeans(lmi, list(pairwise ~ g:time), adjust='none', data=y_mi)[[2]][6])
    theta_extreme_mi_j[j] <- as.numeric(-lsm[2])
    vartotal_extreme_mi_j[j] <- as.numeric(lsm[3]^2)
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
})

stopCluster(cl)

round(rowMeans(result),2)
round(sqrt(apply(result,1,var)),2)
