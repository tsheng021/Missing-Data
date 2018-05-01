library(MASS)
library(nlme)
library(reshape2)
library(parallel)
library(emmeans)

rm(list = ls())


mu_dj <- c(0, 1.3, 2.3, 3.2, 4)
mu_df <- c(0, 1.3, 1.15, 1.6, 2)
mu_pj <- c(0, 1.3, 2.3, 3.2, 4)
mu_pf <- c(0, 1.3, 1.15, 1.6, 2)

sd <- c(2.0, 1.8, 2.0, 2.1, 2.2)
corr <- matrix(c(1.0, 0.6, 0.3, 0.2, 0.1,
                 0.6, 1.0, 0.7, 0.5, 0.2,
                 0.3, 0.7, 1.0, 0.6, 0.4,
                 0.2, 0.5, 0.6, 1.0, 0.5,
                 0.1, 0.2, 0.4, 0.5, 1.0), nrow=5)
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
  y <- mvrnorm(N, mu_dj, vv)
  y <- data.frame(id=1:N, t0=y[,1], t1=y[,2], t2=y[,3], t3=y[,4], t4=y[,5], g=1,
                  stringsAsFactors=FALSE)
  ind <- rep(1, N)
  e <- rep(4, N)
  
  for (j in 2:4){
    cand <- which(ind==1)
    mis <- (j+1):5
    obs <- 1:j
    px <- as.matrix(cbind(1,y[cand,(j+1):(j+2)])) %*% t(t(psi1))
    pt <- 1/(1+exp(-px))
    ind[cand] <- rbinom(length(cand), 1, pt) + ind[cand]
    while (sum(ind[cand]==2)==0) ind[cand] <- rbinom(length(cand), 1, pt) + ind[cand]
    ind_df <- which(ind[cand]==2)
    e[cand][ind_df] <- j
    y[cand,][ind_df, mis+1] <- rnorm(length(ind_df),
                               mean=mu_df[mis]+vv[mis,obs]%*%solve(vv[obs,obs])%*%t(sweep(y[cand,][ind_df,obs+1], 2, mu_dj[obs+1])),
                               sd=sqrt(vv[mis,mis]-vv[mis,obs]%*%solve(vv[obs,obs])%*%vv[obs,mis]))
    
    px <- as.matrix(cbind(1,y[cand,][ind_df,(j+1):(j+2)])) %*% t(t(psi2))
    pt <- 1/(1+exp(-px))
    ind[cand][ind_df] <- ind[cand][ind_df] + rbinom(length(ind_df), 1, pt)
    while (sum(ind[cand]==3)==0|sum(ind[cand]==2)==0) {
      ind[ind[cand]==3] <- ind[ind[cand]==3] - 1
      ind[cand][ind_df] <- ind[cand][ind_df] + rbinom(length(ind_df), 1, pt)
    }
  }
  
  
  # Placebo arm
  z <- mvrnorm(N, mu_pj, vv)
  z <- data.frame(id=1:N+N, t0=z[,1], t1=z[,2], t2=z[,3], t3=z[,4], t4=z[,5], g=0,
                  stringsAsFactors=FALSE)
  indp <- rep(1, N)
  ep <- rep(4, N)
  
  for (j in 2:4){
    cand <- which(indp==1)
    mis <- (j+1):5
    obs <- 1:j
    px <- as.matrix(cbind(1,z[cand,(j+1):(j+2)])) %*% t(t(psi1p))
    pt <- 1/(1+exp(-px))
    indp[cand] <- rbinom(length(cand), 1, pt) + indp[cand]
    while (sum(indp[cand]==2)==0) indp[cand] <- rbinom(length(cand), 1, pt) + indp[cand]
    ind_df <- which(indp[cand]==2)
    ep[cand][ind_df] <- j
    z[cand,][ind_df, mis+1] <- rnorm(length(ind_df),
                                     mean=mu_pf[mis]+vv[mis,obs]%*%solve(vv[obs,obs])%*%t(sweep(z[cand,][ind_df,obs+1], 2, mu_pj[obs+1])),
                                     sd=sqrt(vv[mis,mis]-vv[mis,obs]%*%solve(vv[obs,obs])%*%vv[obs,mis]))
    
    px <- as.matrix(cbind(1,z[cand,][ind_df,(j+1):(j+2)])) %*% t(t(psi2p))
    pt <- 1/(1+exp(-px))
    indp[cand][ind_df] <- indp[cand][ind_df] + rbinom(length(ind_df), 1, pt)
    while (sum(indp[cand]==3)==0|sum(indp[cand]==2)==0) {
      indp[indp[cand]==3] <- indp[indp[cand]==3] - 1
      indp[cand][ind_df] <- indp[cand][ind_df] + rbinom(length(ind_df), 1, pt)
    }
  }
  
  tot_f <- rbind(y, z)
  tot_f$s <- c(ind, indp)
  tot_f$e <- c(e, ep)
  tot_long <- melt(tot_f,id.var=c('id','t0','g','s','e'))
  colnames(tot_long) <- c('id', 'bl', 'g', 's', 'e','time', 'val')

  
  # Create coding variable
  tot_long$s1g1 <- 0
  tot_long$s1g1[which(tot_long$s==1 & tot_long$g==1)] <- 1
  
  tot_long$s23g0 <- 0
  tot_long$s23g0[which(tot_long$s>1 & tot_long$g==0)] <- 1
  
  tot_long$s23g1 <- 0
  tot_long$s23g1[which(tot_long$s>1 & tot_long$g==1)] <- 1

  tot_long$j2s23e2g0 <- 0
  tot_long$j2s23e2g0[which(tot_long$s>1 & 
                          tot_long$g==0 & 
                          tot_long$e==2 &
                          tot_long$time=='t2')] <- 1
  
  tot_long$j2s1g1 <- 0
  tot_long$j2s1g1[which(tot_long$s==1 & 
                             tot_long$g==1 & 
                             tot_long$time=='t2')] <- 1
  
  tot_long$j2s23e2g1 <- 0
  tot_long$j2s23e2g1[which(tot_long$s>1 & 
                             tot_long$g==1 & 
                             tot_long$e==2 &
                             tot_long$time=='t2')] <- 1
  
  tot_long$j3s23e2g0 <- 0
  tot_long$j3s23e2g0[which(tot_long$s>1 & 
                             tot_long$g==0 & 
                             tot_long$e==2 &
                             tot_long$time=='t3')] <- 1
  
  tot_long$j3s23e3g0 <- 0
  tot_long$j3s23e3g0[which(tot_long$s>1 & 
                             tot_long$g==0 & 
                             tot_long$e==3 &
                             tot_long$time=='t3')] <- 1
  
  tot_long$j3s1g1 <- 0
  tot_long$j3s1g1[which(tot_long$s==1 & 
                          tot_long$g==1 & 
                          tot_long$time=='t3')] <- 1
  
  tot_long$j3s23e2g1 <- 0
  tot_long$j3s23e2g1[which(tot_long$s>1 & 
                             tot_long$g==1 & 
                             tot_long$e==2 &
                             tot_long$time=='t3')] <- 1
  
  tot_long$j3s23e3g1 <- 0
  tot_long$j3s23e3g1[which(tot_long$s>1 & 
                             tot_long$g==1 & 
                             tot_long$e==3 &
                             tot_long$time=='t3')] <- 1
  
  tot_long$j4s23e2g0 <- 0
  tot_long$j4s23e2g0[which(tot_long$s>1 & 
                             tot_long$g==0 & 
                             tot_long$e==2 &
                             tot_long$time=='t4')] <- 1
  
  tot_long$j4s23e3g0 <- 0
  tot_long$j4s23e3g0[which(tot_long$s>1 & 
                             tot_long$g==0 & 
                             tot_long$e==3 &
                             tot_long$time=='t4')] <- 1
  
  tot_long$j4s23e4g0 <- 0
  tot_long$j4s23e4g0[which(tot_long$s>1 & 
                             tot_long$g==0 & 
                             tot_long$e==4 &
                             tot_long$time=='t4')] <- 1
  
  tot_long$j4s1g1 <- 0
  tot_long$j4s1g1[which(tot_long$s==1 & 
                          tot_long$g==1 & 
                          tot_long$time=='t4')] <- 1
  
  tot_long$j4s23e2g1 <- 0
  tot_long$j4s23e2g1[which(tot_long$s>1 & 
                             tot_long$g==1 & 
                             tot_long$e==2 &
                             tot_long$time=='t4')] <- 1
  
  tot_long$j4s23e3g1 <- 0
  tot_long$j4s23e3g1[which(tot_long$s>1 & 
                             tot_long$g==1 & 
                             tot_long$e==3 &
                             tot_long$time=='t4')] <- 1
  
  tot_long$j4s23e4g1 <- 0
  tot_long$j4s23e4g1[which(tot_long$s>1 & 
                             tot_long$g==1 & 
                             tot_long$e==4 &
                             tot_long$time=='t4')] <- 1
  
  # MMRM
  m1 <- lme(val~bl+time+bl:time+s1g1+s23g0+s23g1+
              j2s23e2g0+j2s1g1+j2s23e2g1+
              j3s23e2g0+j3s23e3g0+j3s1g1+j3s23e2g1+j3s23e3g1+
              j4s23e2g0+j4s23e3g0+j4s23e4g0+j4s1g1+j4s23e2g1+j4s23e3g1+j4s23e4g1,
            random=~(-1)+time|id,
            method='REML',
            data=tot_long,
            control=lmeControl(returnObject=TRUE))
  
  m1_coef <- summary(m1)$coef$fixed[c(1,3:23)]
  vv_m1 <- vcov(m1)
  sig_b <- vcov(m1)[c(1,3:5),c(1,3:5)]
  sig_g <- vcov(m1)[6:23,6:23]
  cov_bg <- vcov(m1)[6:23, c(1,3:5)]
  
  m2 <- length(which(indp==2 & e==2))
  m3 <- length(which(indp==2 & e==3))
  m4 <- length(which(indp==2 & e==4))
  q2 <- length(which(ind==2 & e==2))
  q3 <- length(which(ind==2 & e==3))
  q4 <- length(which(ind==2 & e==4))
  
  
  sig_o <- matrix(c(m2*(N-m2)/N^3, -m2*m3/N^3,     -m2*m4/N^3,       0,               0,              0, 
                   -m2*m3/N^3,      m3*(N-m3)/N^3, -m3*m4/N^3,       0,               0,              0,
                   -m2*m4/N^3,     -m3*m4/N^3,      m4*(N-m4)/N^3,   0,               0,              0,
                    0,              0,              0,               q2*(N-q2)/N^3,  -q2*q3/N^3,     -q2*q4/N^3,
                    0,              0,              0,              -q2*q3/N^3,       q3*(N-q3)/N^3, -q3*q4/N^3,
                    0,              0,              0,              -q2*q4/N^3,      -q3*q4/N^3,      q4*(N-q4)/N^3
                   ), nrow=6, byrow=T)
  
  ws2 <- length(which(indp>1 & ep==2))/N
  ws3 <- length(which(indp>1 & ep==3))/N
  ws4 <- length(which(indp>1 & ep==4))/N
  wm2 <- length(which(indp==3 & ep==2))/length(which(indp>1 & ep==2))
  wm3 <- length(which(indp==3 & ep==3))/length(which(indp>1 & ep==3))
  wm4 <- length(which(indp==3 & ep==4))/length(which(indp>1 & ep==4))
  
  ps2 <- length(which(ind>1 & e==2))/N
  ps3 <- length(which(ind>1 & e==3))/N
  ps4 <- length(which(ind>1 & e==4))/N
  pm2 <- length(which(ind==3 & e==2))/length(which(ind>1 & e==2))
  pm3 <- length(which(ind==3 & e==3))/length(which(ind>1 & e==3))
  pm4 <- length(which(ind==3 & e==4))/length(which(ind>1 & e==4))
  
  
  ws <- length(which(indp>1))/N
  ps <- length(which(ind>1))/N
  po <- length(which(ind==1))/N
  
  U_marg <- matrix(c(1, 0, 0, 0, 0, 0, 0, 0,
                     1, 0, 1, 0, 0, 0, 0, 0,
                     1, 0, 0, 1, 0, 0, 0, 0,
                     1, 0, 0, 0, 1, 0, 0, 0, 
                     1, 1, 0, 0, 0, 0, 0, 0,
                     1, 1, 1, 0, 0, 1, 0, 0,
                     1, 1, 0, 1, 0, 0, 1, 0,
                     1, 1, 0, 0, 1, 0, 0, 1), nrow=8, byrow=T)
  U_c1 <- matrix(c(1, 0, 0, 0,
                   1, 1, 0, 0,
                   1, 0, 1, 0,
                   1, 0, 0, 1,
                   1, 0, 0, 0,
                   1, 1, 0, 0,
                   1, 0, 1, 0,
                   1, 0, 0, 1), nrow=8, byrow=T)
  u20 <- c(ws2, 0, 0)
  u30 <- c(ws2, ws3, 0, 0, 0)
  u40 <- c(ws2, ws3, ws4, 0, 0, 0, 0)
  u21 <- c(0, po, ps2)
  u31 <- c(0, 0, po, ps2, ps3)
  u41 <- c(0, 0, 0, po, ps2, ps3, ps4)
  
  U_c2 <-  matrix(c(0,    ws, 0,  rep(0, 3), rep(0, 5), rep(0, 7),
                    0,    ws, 0,  u20,       rep(0, 5), rep(0, 7),
                    0,    ws, 0,  rep(0, 3), u30,       rep(0, 7),
                    0,    ws, 0,  rep(0, 3), rep(0, 5), u40, 
                    po,   0,  ps, rep(0, 3), rep(0, 5), rep(0, 7),
                    po,   0,  ps, u21,       rep(0, 5), rep(0, 7),
                    po,   0,  ps, rep(0, 3), u31,       rep(0, 7),
                    po,   0,  ps, rep(0, 3), rep(0, 5), u41 
                    ), nrow=8, byrow=T)
  beta <- matrix(m1_coef[c(1:4)],ncol=1)
  gamma <- matrix(m1_coef[5:22], ncol=1)
  a1 <- solve(t(U_marg)%*%U_marg)%*%t(U_marg)%*%(U_c1%*%beta+U_c2%*%gamma)
  
  o <- matrix(c(m2, m3, m4, q2, q3, q4),
              ncol=1)
  X <- U_c2%*%gamma%*%t(o)%*%ginv(o%*%t(o))
  
  v_a1 <- c(0,1,0,0,0,0,0,1)%*%solve(t(U_marg)%*%U_marg)%*%t(U_marg)%*%
    (X%*%sig_o%*%t(X)+
       U_c1%*%sig_b%*%t(U_c1)+
       U_c2%*%sig_g%*%t(U_c2)+
       U_c1%*%t(cov_bg)%*%t(U_c2)+
       U_c2%*%cov_bg%*%t(U_c1))%*%
    t(solve(t(U_marg)%*%U_marg)%*%t(U_marg))%*%t(t(c(0,1,0,0,0,0,0,1)))
  
  
  
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

