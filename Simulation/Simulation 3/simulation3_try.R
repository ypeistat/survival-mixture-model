rm(list = ls())
source("simulated_data.R")
set.seed(1)
n = 600
K = 2
p = 2
mmax = 10
prop = c(1/3,2/3)
beta = c(1,1)
lambda01 = 1
lambda02 = 0.5
c = 20
data = simulated_data(beta, lambda01, lambda02, c)
###########
lambda <- sqrt(log(n)/n)*0.8
X <- data[,-1:-2]
X <- as.matrix(X)
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

require(survival)
beta0 = lapply(1:mmax, function(k){
  idx = which(class==k)
  result <- coxph(Surv(time, status) ~ X1 + X2, data = data[idx,]) # X1 + X2 need to modify
  return(as.numeric(result$coefficients))
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
# Risk set function
risk.set <- function(t) which(Time >= t)
rs <- apply(as.matrix(Time[as.logical(Delta)]), 1, risk.set)
n.obs <- sum(Delta)
##### 
while(!convergence & iter<iter.max){
  
  ##E-step
  h0 = NULL
  H0 = NULL
  EE = NULL
  PLL = NULL
  allUEXB = NULL
  for (k in 1:mmax){
    Uk = U0[,k]
    XB = X%*%beta0[[k]]  
    EXB = exp(XB)
    EE = cbind(EE, EXB)
    UEXB = Uk*EXB
    allUEXB = cbind(allUEXB, UEXB)
    aa = sapply(Time, function(x)sum(UEXB[Time>=x]))
    bb = 1/aa
    bb[bb == Inf] = 0
    h00 = Uk*bb*Delta
    H00 = sapply(Time, function(x)sum(h00[Time<=x]))
    h0 = cbind(h0, h00)
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
  ####### M step
  if(!convergence){
    U0 = do.call(rbind, lapply(1:n,function(i){
      x = CC[i,]
      return(x / sum(x))
    }))
    prop0 = pmax(0,(colMeans(U0) - lambda) / (1 - mmax*lambda))
    ind <- (prop0>0)
    mmax <- sum(ind)
    ## some problem
    prop0 <- prop0[ind]
    prop0 <- prop0/sum(prop0)
    pdf_est <- pdf_est[,ind]
    CC <- pdf_est %*% diag(as.vector(prop0), mmax)
    U0 = do.call(rbind, lapply(1:n,function(i){
      x = CC[i,]
      return(x / sum(x))
    }))

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
      result <- optim(c(0,0),log.parlik, control = list(maxit = 1000))$par
      return(as.numeric(result))
    }) 
    
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
    class = apply(U0, 1, function(x){
      which.max(x)[1]
    })
    finalmixloglik = mixloglik[length(mixloglik)]
    lkhd = 2*finalmixloglik - (mmax-1 + mmax*p)*log(n)
    print(iter)
  }
}