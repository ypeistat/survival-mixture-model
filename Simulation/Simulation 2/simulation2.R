rm(list = ls())
t1 = Sys.time()
source("fcts.R")
source("cox_EM_nopenalty.R")
source("gendata.R")
source("cox_BIC.R")
set.seed(1)
n = 900
K0 = 2
p = 2
prop = c(2/3,1/3)
beta1 = c(-3,-2)
beta2 = c(1,1)
beta0 = rbind(beta1,beta2)
A = 100
kA_BIC <- rep(0,A)
for (i in 1:A) {
  tryCatch({
    data = gendata(prop,n,K0,p,beta0)
    kA_BIC[i] = cox_BIC(data)
    print(i)
  }, error=function(e){})
}
t2 = Sys.time()
(time = t2 - t1)


hist(kA_BIC[1:100], main="Histogram of estimated numbers of components ",xlim = c(1,10), ylim = c(0,100),
     col = 1:10,xlab = "components")

epan <- function(x){
  return(0.75*(1 - x^2)*((sign(1-x^2)+1)/2))
}
bw = n^(-1/5)
hazard.hat <- function(Time,bw,n,h00){
  h.vec <- rep(0,n)
  for(j in 1:n){
    h.vec[j] <- 1/bw*sum(epan((Time - Time[j])/bw)*h00[j])
  }
  return(h.vec)
}








