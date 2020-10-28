rm(list = ls())
library(glmnet)
library(ggplot2)
library(survminer)
library(PFMC)
source("gendata.R")
source("cox_est.R")
source("cox_EM.R")
set.seed(1)
n = 600
K = 2
p = 2
mmax = 10
prop = c(1/3,2/3)
beta1 = c(-3,-2)
beta2 = c(1,2)
beta0 = rbind(beta1,beta2)
A = 1
kA <- rep(0,A)
kA_IC <- matrix(0, A, 2)
propA = matrix(0, A, mmax)
betaA = list()
lambdaA = rep(0,A)
for (i in 1:A) {
  tryCatch({
  data = gendata(prop,n,K,p,beta0)
  res = cox_est(data,mmax)
  kA[i] = res$k
  k = length(res$pi)
  propA[i,1:k] = res$pi
  betaA[[i]] = res$beta
  lambdaA[i] = res$lambda
  print(i)
  }, error=function(e){})
}

# res <- pmix_cox(Time = data$time,
#                 Delta = data$status,
#                 X = data[,-1:-2],
#                 K = mmax,
#                 tparm = lambda[i])
### No. of components

hist(kA[kA!=0],main="Histogram of estimated numbers of components",
     xlim = c(1,10),ylim = c(0,100),col = c("lightblue", "orange", "gray", "green"),xlab = "components")


### mixing probability
id1 = which(kA==2)
pi0 = propA[id1,]
pi_hat = t(apply(pi0,1,sort)) 
colMeans(pi_hat[,9:10])
sd(pi_hat[,9])
sd(pi_hat[,10])
boxplot(pi_hat[,9:10],col=c("lightblue","orange"),ylim = c(0,1),main="Estimates of mixing probability")



# colMeans(pi_hat[,4:5])
# sd(pi_hat[,4])
# sd(pi_hat[,5])
# boxplot(pi_hat[,4:5],col=c("lightblue","orange"),ylim = c(0,1),main="Estimates of mixing probability")

### beta
#source("Dist_L.R")
beta10 = beta1
beta20 = beta2

beta_est = betaA[id1]
beta1_hat = matrix(0,nrow = length(beta_est),ncol=2)
beta2_hat = matrix(0,nrow = length(beta_est),ncol=2)
#delta1 = rep(0,length(beta_est))
#delta2 = rep(0,length(beta_est))

for (r in 1:length(beta_est)){
  temp1 = unlist(beta_est[[r]][1])
  temp2 = unlist(beta_est[[r]][2])
  # if(norm(temp1,type="2") < norm(temp2,type="2")){
  if(abs(temp1[1]) < abs(temp2[1])){
    tmp1 <- temp2
    tmp2 <- temp1
  }else{
    tmp1 <- temp1
    tmp2 <- temp2
  } 
  beta1_hat[r,] = tmp1
  beta2_hat[r,] = tmp2
}
  #delta1[r] = Dist_L(beta1_hat[r,], beta10)
  #delta2[r] = Dist_L(beta2_hat[r,], beta20)

#res_delta = cbind(delta1,delta2)
# colnames(res_delta)=c("Subgroup 1", "Subgroup2")
# boxplot(res_delta,col=c("lightblue","orange"),ylim = c(0,1))


res_beta1 = beta1_hat
colnames(res_beta1)=c("beta1", "beta2")
boxplot(res_beta1,col=c("lightblue","orange"),main="Subgroup 1",ylim=c(-5,2))

res_beta2 = beta2_hat
colnames(res_beta2)=c("beta1", "beta2")
boxplot(res_beta2,col=c("lightblue","orange"),main="Subgroup 2",ylim=c(0,5))


(bias1 <- colMeans(beta1_hat)- beta1)
(bias2 <- colMeans(beta2_hat) - beta2)
sd(beta1_hat[,1])
sd(beta1_hat[,2])
sd(beta2_hat[,1])
sd(beta2_hat[,2])
  






