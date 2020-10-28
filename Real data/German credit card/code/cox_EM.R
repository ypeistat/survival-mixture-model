cox_EM <- function(data, lambda, mmax){
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
  
  
  # require(survival)
  # beta0 = lapply(1:mmax, function(k){
  #   idx = which(class==k)
  #   result <- coxph(Surv(time, status) ~ X1 + X2 + X3 + X4, data = data[idx,]) # X1 + X2 need to modify
  #   return(as.numeric(result$coefficients))
  # })
  
  
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
      result <- optim(c(0,0,0,0),log.parlik, control = list(maxit = 1000))$par
      return(as.numeric(result))
    }) 
    
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
      #bb[bb == Inf] = 0
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
    ### E step 
    U0 = do.call(rbind, lapply(1:n,function(i){
      x = CC[i,]
      return(x / sum(x))
    }))
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
  ####
  result = list(U = U0,
                k = mmax,
                CC = CC,
                pi = prop0,
                beta = beta0,
                BIC = BIC,
                class = class,
                ploglik = ploglik,
                mixloglik = mixloglik,
                iter = iter,
                convergence = iter < iter.max)
  return(result)
}

#plot(ploglikpenalty,type = 'o',col='blue',xlab="iter", main = "Penalized log-likelihood function")



 