library(MASS)
library(clime)
library(tobitnet)

source('script/TobitEstimate.R')
source('script/ghFunction.R')


InferenceTheta <- function(x,y,c_q,alpha,B){
  n <- nrow(x)
  p <- ncol(x)
  d <- as.numeric(y > c_q) 
  thetahat <- theta_tobit(x, y, c_q)
  delta <- thetahat[1:p]
  gamma <- thetahat[(p+1)] 
  
  J_11 <- matrix(0,p,p)
  J_12 <- matrix(0, p, 1)
  J_22 <- 0
  for(i in 1:n){
    x_I <- matrix(x[i, ], ncol = 1) 
    xI_delta <- as.numeric(t(x_I)%*%delta)
    s <- gamma*c_q-xI_delta
    gh <- gh_stable(s)
    g_s <- gh$g
    h_s <- gh$h
    J_11 <- J_11 + (d[i]+(1-d[i])*h_s)*(x_I%*%t(x_I)) 
    J_12 <- J_12 + ( -d[i]*y[i] - (1-d[i])*c_q*h_s ) * x_I 
    J_22 <- J_22 + d[i]*( y[i]^2 + gamma^(-2) ) + (1-d[i])*( c_q^2 * h_s )
  }
  J_11 <- J_11/n
  J_12 <- J_12/n
  J_22 <- J_22/n
  J_21 <- t(J_12)
  J <- rbind(cbind(J_11, J_12), c(J_21, J_22)) 
  
  C2 <- 0.4  #The range is from 0.4 to 1.2
  J_inverse_hat <- clime(J,lambda=C2* sqrt(log(p + 1) / n),sigma=TRUE,linsolver="simplex")
  W_hat <- J_inverse_hat$Omega[[1]]
  W_hat <- (W_hat + t(W_hat))/2
  
  first_order_1 <- matrix(0,p,1)
  first_order_2 <- 0
  for(i in 1:n){
    x_I <- matrix(x[i, ], ncol = 1)
    xI_delta <- as.numeric(t(x_I)%*%delta)
    s <- gamma*c_q-xI_delta
    gh <- gh_stable(s)
    g_s <- gh$g
    first_order_1 <- first_order_1 + (-d[i]*x_I*(gamma*y[i]-xI_delta))+(1-d[i])*x_I*g_s
    first_order_2 <- first_order_2 + d[i]*( -1/gamma + y[i]*(gamma*y[i]-xI_delta) ) - (1-d[i])*c_q*g_s
  }
  first_order_1 <- first_order_1/n
  first_order_2 <- first_order_2/n
  first_order <- -c(first_order_1, first_order_2)
  debias <- as.numeric(thetahat)+ as.numeric(W_hat%*%first_order)
  
  S_hat_sqrt <- sqrt(diag(W_hat))
  S_hat_sqrt_inverse <- matrix(0,(p+1),(p+1))
  diag(S_hat_sqrt_inverse) <- 1/S_hat_sqrt  
  bootstraps <- sapply(1:B,function(xxx,x,y,delta,gamma,W_hat,S_hat_sqrt_inverse){
    n <- nrow(x)
    p <- ncol(x)
    e_hat = rnorm(n)
    
    first_order_1e <- matrix(0,p,1)
    first_order_2e <- 0
    for(i in 1:n){
      x_I <- matrix(x[i, ], ncol = 1)
      xI_delta <- as.numeric(t(x_I)%*%delta)
      s <- gamma*c_q-xI_delta
      gh <- gh_stable(s)
      g_s <- gh$g
      first_order_1e <- first_order_1e + e_hat[i]*((-d[i]*x_I*(gamma*y[i]-xI_delta))+(1-d[i])*x_I*g_s)
      first_order_2e <- first_order_2e + e_hat[i]*(d[i]*( -1/gamma + y[i]*(gamma*y[i]-xI_delta) ) - (1-d[i])*c_q*g_s)
    }
    first_order_1e <- first_order_1e/n
    first_order_2e <- first_order_2e/n
    first_order_e <- -c(first_order_1e, first_order_2e)
    boots <- sqrt(n)*S_hat_sqrt_inverse%*%W_hat%*%first_order_e
    return(max(abs(boots)))
  },x,y,delta,gamma,W_hat,S_hat_sqrt_inverse)
  Q <- as.numeric(quantile(bootstraps, probs = 1 - alpha, names = FALSE)) / sqrt(n)
  lowerbound <- debias - as.numeric(S_hat_sqrt)*Q
  upperbound <- debias + as.numeric(S_hat_sqrt)*Q
  return(list(lowerbound = lowerbound, upperbound = upperbound))
}

InferenceBeta <- function(x,y,c_q,alpha,B){
  n <- nrow(x)
  p <- ncol(x)
  d <- as.numeric(y > c_q) 
  thetahat <- theta_tobit(x, y, c_q)
  delta <- thetahat[1:p]
  gamma <- thetahat[(p+1)] 
  
  J_11 <- matrix(0,p,p)
  J_12 <- matrix(0, p, 1)
  J_22 <- 0
  for(i in 1:n){
    x_I <- matrix(x[i, ], ncol = 1) 
    xI_delta <- as.numeric(t(x_I)%*%delta)
    s <- gamma*c_q-xI_delta
    gh <- gh_stable(s)
    g_s <- gh$g
    h_s <- gh$h
    J_11 <- J_11 + (d[i]+(1-d[i])*h_s)*(x_I%*%t(x_I)) 
    J_12 <- J_12 + ( - d[i]*y[i] - (1-d[i]) * c_q * h_s ) * x_I 
    J_22 <- J_22 + d[i]*((y[i])^2+gamma^(-2)) + (1-d[i]) * ( c_q^2 * h_s )
  }
  J_11 <- J_11/n
  J_12 <- J_12/n
  J_22 <- J_22/n
  J_21 <- t(J_12)
  J <- rbind(cbind(J_11, J_12), c(J_21, J_22)) 
  
  C2 <- 0.4  #The range is from 0.4 to 1.2
  J_inverse_hat <- clime(J,lambda=C2*sqrt(log(p + 1) / n),sigma=TRUE,linsolver="simplex")
  W_hat <- J_inverse_hat$Omega[[1]]
  W_hat <- (W_hat + t(W_hat))/2
  
  first_order_1=matrix(0,p,1)
  first_order_2=0
  for(i in 1:n){
    x_I <- matrix(x[i, ], ncol = 1)
    xI_delta <- as.numeric(t(x_I)%*%delta)
    s <- gamma*c_q-xI_delta
    gh <- gh_stable(s)
    g_s <- gh$g
    first_order_1 = first_order_1 + (-d[i]*x_I*(gamma*y[i]-xI_delta))+(1-d[i])*x_I*g_s
    first_order_2 = first_order_2 + d[i] * ( -1/gamma + y[i] * (gamma*y[i] - xI_delta) ) -(1-d[i]) * c_q * g_s
  }
  first_order_1 <- first_order_1/n
  first_order_2 <- first_order_2/n
  first_order <- -c(first_order_1, first_order_2)
  
  debias <- thetahat+W_hat%*%first_order
  delta_t <- as.numeric(debias[1:p])
  gamma_t <- as.numeric(debias[p+1])
  betadebias <- delta_t / gamma_t
  
  S_hat_sqrt <- sqrt(diag(W_hat))
  S_hat_sqrt_inverse <- matrix(0,(p+1),(p+1))
  diag(S_hat_sqrt_inverse) <- 1/S_hat_sqrt  
  
  bootstraps <- sapply(1:B,function(xxx,x,y,delta,gamma,W_hat,S_hat_sqrt_inverse){
    n <- nrow(x)
    p <- ncol(x)
    e_hat <- rnorm(n)
    
    first_order_1e <- matrix(0,p,1)
    first_order_2e <- 0
    for(i in 1:n){
      x_I <- matrix(x[i, ], ncol = 1)
      xI_delta <- as.numeric(t(x_I)%*%delta)
      s <- gamma*c_q-xI_delta
      gh <- gh_stable(s)
      g_s <- gh$g
      first_order_1e <- first_order_1e + e_hat[i]*((-d[i]*x_I*(gamma*y[i]-xI_delta))+(1-d[i])*x_I*g_s)
      first_order_2e <- first_order_2e + e_hat[i]*( d[i] * ( -1/gamma + y[i] * (gamma*y[i] - xI_delta) ) - (1-d[i]) * c_q * g_s )
    }
    first_order_1e <- first_order_1e/n
    first_order_2e <- first_order_2e/n
    first_order_e <- -c(first_order_1e, first_order_2e)
    
    boots <- sqrt(n)*W_hat%*%first_order_e
    bootdelta <- boots[1:p]; bootgamma <- boots[(p+1)]
    Shatsqrtinverse <- S_hat_sqrt_inverse[1:p,1:p]
    bootbeta <- Shatsqrtinverse%*%(bootdelta/gamma-delta*bootgamma/(gamma^2))
    return(max(abs(bootbeta[1:p])))
  },x,y,delta,gamma,W_hat,S_hat_sqrt_inverse)
  Q <- as.numeric(quantile(bootstraps, probs = 1 - alpha, names = FALSE)) / sqrt(n)
  SS <- S_hat_sqrt[1:p]
  lowerbound <- betadebias - SS*Q
  upperbound <- betadebias + SS*Q
  return(list(lowerbound = lowerbound, upperbound = upperbound))
}






