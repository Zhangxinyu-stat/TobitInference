library(tobitnet)

theta_tobit <- function(x,y,c_q) {
  n <- nrow(x)
  p <- ncol(x)
  C1 <- 1.6  #The range is from 0.8 to 2
  lambda <- C1 * sqrt(log(p) / n)  
  tb1 <- tobitscad(x, y, c = c_q, lambda = lambda)
  sigma_hat <- tb1$sigma   
  gamma_hat <- 1/sigma_hat
  delta_hat <- as.numeric(tb1$beta) * gamma_hat   
  theta_hat <- c(delta_hat, gamma_hat)
  return(theta_hat)
}