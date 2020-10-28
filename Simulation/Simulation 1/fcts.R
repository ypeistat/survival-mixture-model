bh.smooth <- function(Time,bw,n,h00){
  h.vec <- rep(0,n)
  for(i in 1:n){
    h.vec[i] <- 1/bw*sum(phi((Time - Time[i])/bw)*h00)
  }
  return(h.vec)
}

epan <- function(x){
  return(0.75*(1 - x^2)*((sign(1-x^2)+1)/2))
}

phi <- function(x){
  exp(-0.5*x^2)/sqrt(2*pi)
}