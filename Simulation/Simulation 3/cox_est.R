cox_est <- function(data,mmax){
  n <- dim(data)[1]
  lambda <- seq(sqrt(log(n)/n)*0.1, sqrt(log(n)/n)*0.5,length.out = 5)
  #lambda <- sqrt(log(n)/n)*1.0 ## mmax = 5
  #lambda <- sqrt(log(n)/n)*0.31 ## mmax = 10
  lkhd_max <- -10^6
  for (i in 1:length(lambda)){
     res <- cox_EM(data,lambda[i],mmax)
     if(lkhd_max < res$lkhd){
       k_max <- res$k
       prop_max <- res$pi
       beta_max <- res$beta
       lkhd_max <- res$lkhd
       lambda_max <- lambda[i]
     }
  }
  result <- list(pi = prop_max,
                 k = k_max,
                 beta = beta_max,
                 lkhd  = lkhd_max,
                 lambda = lambda_max)
  return(result)
}

#plot(ploglikpenalty,type = 'o',col='blue',xlab="iter", main = "Penalized log-likelihood function")