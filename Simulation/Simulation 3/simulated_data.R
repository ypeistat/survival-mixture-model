simulated_data <- function(beta, lambda01, lambda02, c){
  X = matrix(0,nrow = n,ncol=p) ###predictors
  X[,1] = rnorm(n)
  rho <- 0.5
  for(j in 2:p){
    X[,j] = rho*X[,j-1] + sqrt(1-rho^2)*rnorm(n)
  } 
  
  lambda1 <- lambda01 * exp(X %*% beta)
  lambda2 <- lambda02 * exp(X %*% beta)
  t1 <- rexp(n*prop[1], rate = lambda1)
  t2 <- rexp(n*prop[2], rate = lambda2)
  time <- c(t1,t2)
  C <- runif(n, 0, c)
  tildeT <- apply(cbind(C, time), 1, 'min')
  delta <- (time < C)
  N <- sum(delta==1) # Number of observed events
  data = as.data.frame(cbind( tildeT,delta,X))
  colnames(data) = c("time", "status","X1","X2")
  return(data)
}

