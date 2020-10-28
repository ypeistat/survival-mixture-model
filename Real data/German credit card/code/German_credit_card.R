################## glmnet select variable cox model
rm(list=ls())
library(knitr)
library(dplyr)
library(tidyr)
library(reshape2)
library(RColorBrewer)
library(GGally)
library(ggplot2)
library(caret)
library(glmnet)
library(boot)
#library(verification)
library(survival)
#library(condSURV)
#library(JM)
library(survminer)
library(survivalROC)

german_credit <- read.table("http://archive.ics.uci.edu/ml/machine-learning-databases/statlog/german/german.data")

colnames(german_credit) <- c("chk_acct", "duration", "credit_his", "purpose", 
                             "amount", "saving_acct", "present_emp", "installment_rate", "sex", "other_debtor", 
                             "present_resid", "property", "age", "other_install", "housing", "n_credits", 
                             "job", "n_people", "telephone", "foreign", "status")
german_credit$status <- german_credit$status - 1

##### split data into training data and test data
set.seed(123456)
in.train <- createDataPartition(as.factor(german_credit$status), p=0.8, list=FALSE)
german_credit.train <- german_credit[in.train,]
german_credit.test <- german_credit[-in.train,]

### transform categorical variable to factor
factor_var <- c(1,3,4,6,7,9,10,12,14,15,17,19,20)
num_var <- c(5,8,11,13,16,18)
train <- german_credit.train
test <- german_credit.test
train[num_var] <- scale(train[num_var])
test[num_var] <- scale(test[num_var])
train[factor_var] <- sapply(train[factor_var] , as.numeric)
test[factor_var] <- sapply(test[factor_var] , as.numeric)

### prepare train data to glmnet
X.train <- as.matrix(train[,c(1,3:20)])
time <- train[,2]
status <- train[,21]
surv.train = cbind(time = time, status = status)
#######
X.test <- as.matrix(test[,c(1,3:20)])
time1 <- test[,2]
status1 <- test[,21]
surv.test = cbind(time = time1, status = status1)
### fit glmnet to select important features
result = glmnet(x=X.train, y= surv.train, family = c("cox"))
plot(result, xvar = "lambda", label=TRUE)
cv.lasso<- cv.glmnet(x=X.train, y=surv.train, family = "cox", alpha = 1, nfolds = 10)
plot(cv.lasso)
cv.lasso$lambda.1se
coef(result, s=cv.lasso$lambda.1se)

####################################
### final data to be used for further analysis
sig_id <- c(1,5)
df <- cbind(surv.train,train[,sig_id])
df_test <- cbind(surv.test,test[,sig_id])
own_check_id <- which(df[,3]!=4)
own_check_id1 <- which(df_test[,3]!=4)
df[own_check_id,3]=1
df_test[own_check_id1,3]=1
df[-own_check_id,3]=0
df_test[-own_check_id1,3]=0
data <- df
colnames(data) = c("time", "status","X1","X2")

################################# survival mixture model
source("fcts.R")
mmax = 10
n = dim(data)[1]
#lambda = seq(sqrt(log(n)/n)/30,sqrt(log(n)/n)/15,length.out = 16)
lambda = sqrt(log(n)/n)/20


#########################
X <- data[,-1:-2]
X <- as.matrix(X)
n <- dim(X)[1]
p <- dim(X)[2]
Time = matrix(data[,1], ncol = 1)
Delta = matrix(data[,2], ncol = 1)
surv = cbind(time = Time, status = Delta)
colnames(surv) = c("time", "status")
#### initialization
u.init = sample(1:mmax, size = n, replace = TRUE)
U0 = matrix(0.001, mmax, mmax)
diag(U0) = 1 - (mmax-1)*0.001
U0 = U0[u.init, ]
#####initial pi                                                                                                        
prop0 = apply(U0, 2, mean)
#### initial class label
class = apply(U0, 1, function(x){
  which.max(x)[1]
})


ploglik = NULL
mixloglik = NULL
ploglikpenalty = NULL
mixloglikpenalty = NULL
lastploglikpenalty = -Inf
lastmixloglikpenalty = -Inf
convergence = FALSE
iter = 0
iter.max = 1000
epsilon = 1e-6
abstol = 1e-4
reltol = 1e-6
##### 
while(!convergence & iter<iter.max){
  # weighted log partial likelihood function
  beta0 = lapply(1:mmax, function(k){
    Uk <- U0[,k]
    log.parlik <- function(beta){
      XB = X%*%beta  
      EXB = exp(XB)
      UEXB = Uk*EXB
      aa = sapply(Time, function(x)sum(UEXB[Time>=x]))
      bb = 1/aa
      bb[bb == Inf] = 0
      PLL = sum((Uk*(XB-log(aa)))[Delta==1&aa>0])
      return(-PLL)
    }
    result <- optim(c(0,0),log.parlik, control = list(maxit = 1000),hessian = TRUE)
    betahat <- result$par
    se <- sqrt(diag(solve(result$hessian)))
    return(as.numeric(cbind(betahat,se)))
  }) 
  
  ##E-step
  h0 = NULL
  H0 = NULL
  EE = NULL
  PLL = NULL
  allUEXB = NULL
  for (k in 1:mmax){
    Uk = U0[,k]
    XB = X%*%beta0[[k]][1:p]  
    EXB = exp(XB)
    EE = cbind(EE, EXB)
    UEXB = Uk*EXB
    allUEXB = cbind(allUEXB, UEXB)
    aa = sapply(Time, function(x)sum(UEXB[Time>=x]))
    bb = 1/aa
    bb[bb == Inf] = 0
    h00 = Uk*bb*Delta
    H00 = sapply(Time, function(x)sum(h00[Time<=x]))
    h0_smooth  = bh.smooth(Time, bw = IQR(Time)*n^(-1/5),n,h00)
    h0 = cbind(h0, h0_smooth)
    H0 = cbind(H0, H00)
    PLL = c(PLL, sum((Uk*(XB-log(aa)))[Delta==1&aa>0]))
  }
  
  pdf_est = (h0 * EE)^(Delta %*% t(rep(1, mmax))) * exp(-H0 * EE)
  CC = pdf_est %*% diag(as.vector(prop0), mmax)
  
  nowploglik = sum(PLL) + sum((U0%*%diag(log(prop0), mmax))[U0 != 0])
  nowmixloglik = sum(log(apply(CC, 1, sum)))
 
  ####
  penalty = function(prop) {
    penalty = n * lambda  * sum(log((epsilon + prop) / epsilon))
    return(penalty)
  }
  
  nowploglikpenalty = nowploglik - penalty(prop0)
  nowmixloglikpenalty = nowmixloglik - penalty(prop0)
  ploglikpenalty = c(ploglikpenalty, nowploglikpenalty)
  mixloglikpenalty = c(mixloglikpenalty, nowmixloglikpenalty)
  ploglik = c(ploglik, nowploglik)
  mixloglik = c(mixloglik, nowmixloglik)
  ############################
  if(is.finite(nowploglikpenalty) & is.finite(lastploglikpenalty) &
     (abs(nowploglikpenalty-lastploglikpenalty) < abstol |
      abs(nowploglikpenalty-lastploglikpenalty) <
      abs(reltol*(nowploglikpenalty + reltol))))
    convergence = TRUE
  
  lastploglikpenalty = nowploglikpenalty
  lastmixloglikpenalty = nowmixloglikpenalty
  iter = iter + 1
  ######
  if(iter==iter.max)
    warning("EM iterations reach iter.max!")
  ####### M step
  if(!convergence){
    U0 = do.call(rbind, lapply(1:n,function(i){
      x = CC[i,]
      return(x / sum(x))
    }))
    prop0 = pmax(0,(colMeans(U0) - lambda) / (1 - mmax*lambda))
    ind <- (prop0>0)
    mmax <- sum(ind)
    ########################
    prop0 <- prop0[ind]
    prop0 <- prop0/sum(prop0)
    pdf_est <- pdf_est[,ind]
    CC <- pdf_est %*% diag(as.vector(prop0), mmax)
    U0 = do.call(rbind, lapply(1:n,function(i){
      x = CC[i,]
      return(x / sum(x))
    }))
  }
  ########## prediction
  class = apply(U0, 1, function(x){
    which.max(x)[1]
  })
  finalmixloglik = mixloglik[length(mixloglik)]
  # cn = 5*log(log(n+mmax))
  BIC = 2*finalmixloglik - (mmax-1 + mmax*p)*log(n)
  print(iter)
  print(mmax)
  print(prop0)
  print(beta0)
  print(class)
}

#################  Evaluate mixture AUC training data 
lp1 <- (X - colMeans(X))%*%beta0[[1]][1:p]
lp2 <- (X - colMeans(X))%*%beta0[[2]][1:p]
df.train <- data
df.train$lp1 <- lp1
df.train$lp2 <- lp2
df.train$lp.mix <- prop0[1]*lp1 + prop0[2]*lp2

fun.survivalROC.train <- function(lp, t) {
  res <- with(df.train,
              survivalROC(Stime        = time,
                          status       = status,
                          marker       = get(lp),
                          predict.time = t,
                          method       = "KM"))       # KM method without smoothing
  
  ## Plot ROCs
  with(res, plot(TP ~ FP, type = "l", main = sprintf("t = %.0f, AUC = %.2f", t, AUC),col=t))
  abline(a = 0, b = 1, lty = 2)
  res
}
layout(matrix(1:6, byrow = T, ncol = 3))
## Model with training data
res.survivalROC.train <- lapply(4:9* 6, function(t) {
  fun.survivalROC.train(lp = "lp.mix", t)
})
################################# no group structure
train.fit = coxph(Surv(time, status) ~ X1 + X2, data = df.train, method = "breslow")
summary(train.fit)
df.train$lp.cph <- predict(train.fit, type = "lp")
res.survivalROC.train <- lapply(4:9 * 6, function(t) {
  fun.survivalROC.train(lp = "lp.cph", t)
})
###################################### AUC test data
test.data <- df_test[,3:4]
newX <- test.data - colMeans(test.data)
lp_test1 <- newX[,1]*beta0[[1]][1] + newX[,2]*beta0[[1]][2]
lp_test2 <- newX[,1]*beta0[[2]][1] + newX[,2]*beta0[[2]][2]
df_test$lp_test1 <- lp_test1
df_test$lp_test2 <- lp_test2
df_test$lp.mix <- prop0[1]*lp_test1 + prop0[2]*lp_test2

fun.survivalROC.test <- function(lp, t) {
  res <- with(df_test,
              survivalROC(Stime        = time,
                          status       = status,
                          marker       = get(lp),
                          predict.time = t,
                          method       = "KM"))       # KM method without smoothing
  ## Plot ROCs
  with(res, plot(TP ~ FP, type = "l", main = sprintf("t = %.0f, AUC = %.2f", t, AUC),col = t))
  abline(a = 0, b = 1, lty = 2)
}
layout(matrix(1:6, byrow = T, ncol = 3))
## Model with test data
res.survivalROC.test <- lapply(4:9* 6, function(t) {
  fun.survivalROC.test(lp = "lp.mix", t)
})
######## no group test data
df_test$lp.cph  <- newX[,1]*train.fit$coefficients[1] + newX[,2]*train.fit$coefficients[2]
res.survivalROC.test <- lapply(4:9* 6, function(t) {
  fun.survivalROC.test(lp = "lp.cph", t)
})


###################################
# df2 <- as.data.frame(cbind(class,data))
# ggplot(df2, aes(class, ..count..)) + 
#   geom_bar(aes(fill = as.factor(status)), position = "dodge")


# model1 <- coxph(Surv(time, status) ~ X2  + X1 + strata(as.factor(class)), data = data)
# ggsurvplot(survfit(model1), data = data, conf.int = FALSE)


################## explaination
colnames(data) = c("time", "status","check_account","X2")
model2 <- coxph(Surv(time, status) ~ X2 + strata(check_account), data = data)
ggsurvplot(survfit(model2), data = data, conf.int = TRUE)
################


##########################prediction
Xnew = df_test[,3:4]
Xnew <- as.matrix(Xnew)
Time_test = matrix(df_test[,1], ncol = 1)
Delta_test = matrix(df_test[,2], ncol = 1)
EE = NULL
for (k in 1:2){
  XB = Xnew%*%beta0[[k]][1:p]  
  EXB = exp(XB)
  EE = cbind(EE, EXB)
}
## assume the baseline hazard function is same, so we can cancel it
pdf_test = EE^(Delta_test %*% t(rep(1, 2))) * exp(EE)
CC_test = pdf_test %*% diag(as.vector(prop0), 2)
U0_est = do.call(rbind, lapply(1:200,function(i){
  x = CC_test[i,]
  return(x / sum(x))
}))

ratio <- U0_est[,2]/U0_est[,1]
g2 <- which(ratio>1)  ## second group
g1 <- which(ratio<1)  ## first group


df_test1 <- cbind(df_test[g1,1:4],rep(1,length(g1)))
colnames(df_test1) <- c("time", "status", "chk_acct", "amount", "class")
                            
df_test2 <- cbind(df_test[g2,1:4],rep(2,length(g2)))
colnames(df_test2) <- c("time", "status", "chk_acct", "amount", "class")
df_test_all <- rbind(df_test1,df_test2)


model_test1 <- coxph(Surv(time, status) ~ strata(class), data = df_test_all)
ggsurvplot(survfit(model_test1), data = df_test_all, conf.int = TRUE)









