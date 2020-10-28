cox_BIC <- function(data){
  bic <- rep(0,10)
  #lkhd <- rep(0,10)
  for (j in 1:10){
    nopenalty.fit = cox_EM_nopenalty(data, K=j)
    bic[j] = nopenalty.fit$BIC
    #lkhd[j] = nopenalty.fit$finalmixloglik
    print(j)
  }
  #print(bic)
  khat_bic <- which.max(bic)
return(khat_bic)
}

# bic_cox <- function(result){
#   iter <- result$iter
#   finalmixloglik = result$mixloglik
#   pi_hat = result$pi
#   Khat = length(pi_hat)
#   beta = result$beta
#   p = length(beta)
#   n = dim(result$U)[1]
#   BIC = 2*finalmixloglik[iter] - (Khat-1+Khat*p)*log(n)
#   return(BIC)
# }


