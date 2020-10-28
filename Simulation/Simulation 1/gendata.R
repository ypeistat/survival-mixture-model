gendata <- function(prop,n,K,p,beta0){
  c = 5 #control censoring proportions: 5%, 5.0; 25%, 1
  X = NULL
  time = NULL
  delta = NULL
  group = NULL
  for (k in 1:K){
    nk = n*prop[k]
    xk = matrix(0,nrow = nk,ncol=p) ###predictors
    xk[,1] = rnorm(nk)
    rho <- 0.5
    for(j in 2:p){
      xk[,j] = rho*xk[,j-1]+sqrt(1-rho^2)*rnorm(nk)
    }
    xk = as.matrix(xk)
    X = rbind(X,xk)
    temp = rexp(nk)
    dtime = as.numeric(0.25*log(0.5*temp*exp(-xk%*%beta0[k,])+1.0))
    ctime = c*runif(nk)
    time0 = apply(cbind(ctime,dtime),1,min)
    time = rbind(time,as.matrix(time0))
    delta0 = as.numeric(dtime <= ctime)
    delta = rbind(delta,as.matrix(delta0))
    group0 = rep(k,nk)
    group = c(group,group0)
  }
  data = as.data.frame(cbind(time,delta,X))
  colnames(data) = c("time", "status","X1","X2")
  return(data)
}


# ## plot Kaplan Meier
# simulated_data <- data
# simulated_data$group <- as.factor(group)
# model <- survfit(Surv(time, status) ~ group, data = simulated_data)
# survminer::ggsurvplot(model, conf.int = FALSE)


