library(MASS)
library(reshape2)

mu <- c(1, 3)
sd <- c(1.8, 2.2)
corr <- matrix(c( 1.0, 0.2,
                  0.2, 1.0), nrow=2)
vv <- diag(sd) %*% corr %*% diag(sd)
N <- 100
S <- 5000
k_range <- c(-2, -1.5, -.5, 0, .5, 1, 1.5, 2)
# k_range <- c(-1, -.8, -.6 -.4, -.2, 0, .2, .4, .6, .8, 1)
# k_range <- c(-3.2, -1.6, -0.8, -0.4, -0.2, 0, .2, .4, .8, 1.6, 3.2)
# k_range <- c(-3.2,  3.2)
curve <- c()
for(k in k_range){
  no_cores <- detectCores()-1
  cl <- makeCluster(no_cores)
  clusterExport(cl,varlist=c('S', 'N', 'vv', 'mu', 'k'))
  clusterSetRNGStream(cl, 1234)
  clusterEvalQ(cl, {
    library(MASS)
    library(nlme)
    library(reshape2)
  })
  
  result <- parSapply(cl, 1:S, function(x){
    out <- c()
    y <- data.frame(id=1:N, mvrnorm(N, mu, vv)) # Simulate complete
    colnames(y) <- c('id', 't1', 't2')
    psi <- c(k, -.2, 0)
    # psi <- c(k, -k, 0)
    px <- as.matrix(cbind(1,y[,2:3])) %*% t(t(psi))
    pt <- 1/(1+exp(-px))
    ind <- rbinom(N, 1, pt)+1
    if (sum(ind==2)==0) ind <- rbinom(N, 1, pt)+1
    ind <- which(ind==2) # Simulate dropout MAR
    # ind <- which(y$id>50) # Simulate dropout MCAR
    y$t2[ind] <- NA
    y_mis_long <- na.omit(melt(y,id.vars = 'id'))
    colnames(y_mis_long) <- c('id', 'time', 'val')
    
    m1 <- lme(val~time, 
              random=~(-1)+time|id, 
              method='REML', 
              data=y_mis_long)
    coeff <- m1$coeff$fixed
    D <- c(1,1)
    mu_hat <- D%*%coeff
    sd_hat <- sqrt(D%*%vcov(m1)%*%D)
    out <- c(out, mu_hat, sd_hat)
    
    # MI
    MI <- 20
    sigma <- getVarCov(m1) + diag(2)*m1$sigma^2
    theta_mar_mi_j <- vartotal_mar_mi_j <- rep(NA, MI)
    for(j in 1:MI){
      z_act <- matrix(c(1,0,
                        1,1),nrow=2,byrow=T)
      mu_act <- z_act %*% coeff
      mu_act_tilde <- mvrnorm(1, mu_act, sigma/(N-length(ind)))
      sigma_tilde <- matrix(rWishart(1, N-length(ind)-1, sigma),nrow=2)/(N-length(ind)-1)
      y$t2[ind] <- rnorm(length(ind), 
                         mu_act_tilde[2]+sigma_tilde[2,1]/sigma_tilde[1,1]*(y$t1[ind]-mu_act_tilde[1]),
                         sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
      
      theta_mar_mi_j[j] <- mean(y$t2)
      vartotal_mar_mi_j[j] <- var(y$t2)/N
    }
    theta_mar_mi <- mean(theta_mar_mi_j)
    W <- mean(vartotal_mar_mi_j)
    B <- var(theta_mar_mi_j)
    se_mar_mi <-  sqrt(W + (1+1/MI)*B)
    out <- c(out, theta_mar_mi, se_mar_mi)
  })
  stopCluster(cl)
  
  curve <- rbind(curve, 
                 c(k, round(rowMeans(result),2)[c(2,4)], 
                   round(sqrt(apply(result,1,var)),2)[c(1,3)]))
}


plot(curve[,4]~curve[,1],type='l')
lines(curve[,2]~curve[,1],col='red')
lines(curve[,3]~curve[,1],col='red')
lines(curve[,5]~curve[,1])
# no_cores <- detectCores()-1
# cl <- makeCluster(no_cores)
# clusterExport(cl,varlist=c('S', 'N', 'vv', 'mu', 'k'))
# clusterSetRNGStream(cl, 1234)
# clusterEvalQ(cl, {
#   library(MASS)
#   library(nlme)
#   library(reshape2)
# })
# 
# result <- parSapply(cl, 1:S, function(x){
#   out <- c()
#   y <- data.frame(id=1:N, mvrnorm(N, mu, vv)) # Simulate complete
#   colnames(y) <- c('id', 't1', 't2')
#   psi <- c(k, -k, 0)
#   px <- as.matrix(cbind(1,y[,2:3])) %*% t(t(psi))
#   pt <- 1/(1+exp(-px))
#   ind <- rbinom(N, 1, pt)+1
#   if (sum(ind==2)==0) ind <- rbinom(N, 1, pt)+1
#   ind <- which(ind==2) # Simulate dropout MAR
#   # ind <- which(y$id>50) # Simulate dropout MCAR
#   y$t2[ind] <- NA
#   y_mis_long <- na.omit(melt(y,id.vars = 'id'))
#   colnames(y_mis_long) <- c('id', 'time', 'val')
#   
#   m1 <- lme(val~time, 
#             random=~(-1)+time|id, 
#             method='REML', 
#             data=y_mis_long)
#   coeff <- m1$coeff$fixed
#   D <- c(1,1)
#   mu_hat <- D%*%coeff
#   sd_hat <- sqrt(D%*%vcov(m1)%*%D)
#   out <- c(out, mu_hat, sd_hat)
#   
#   # MI
#   MI <- 20
#   sigma <- getVarCov(m1) + diag(2)*m1$sigma^2
#   theta_mar_mi_j <- vartotal_mar_mi_j <- rep(NA, MI)
#   for(j in 1:MI){
#     z_act <- matrix(c(1,0,
#                       1,1),nrow=2,byrow=T)
#     mu_act <- z_act %*% coeff
#     mu_act_tilde <- mvrnorm(1, mu_act, sigma/(N-length(ind)))
#     sigma_tilde <- matrix(rWishart(1, N-length(ind)-1, sigma),nrow=2)/(N-length(ind)-1)
#     y$t2[ind] <- rnorm(length(ind), 
#                            mu_act_tilde[2]+sigma_tilde[2,1]/sigma_tilde[1,1]*(y$t1[ind]-mu_act_tilde[1]),
#                        sqrt(sigma_tilde[2,2]-sigma_tilde[1,2]^2/sigma_tilde[1,1]))
# 
#     theta_mar_mi_j[j] <- mean(y$t2)
#     vartotal_mar_mi_j[j] <- var(y$t2)/N
#   }
#   theta_mar_mi <- mean(theta_mar_mi_j)
#   W <- mean(vartotal_mar_mi_j)
#   B <- var(theta_mar_mi_j)
#   se_mar_mi <-  sqrt(W + (1+1/MI)*B)
#   out <- c(out, theta_mar_mi, se_mar_mi)
# })
# stopCluster(cl)
# 
# round(rowMeans(result),2)
# round(sqrt(apply(result,1,var)),2)
