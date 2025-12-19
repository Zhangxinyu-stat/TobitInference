gh_stable <- function(s){
  logPhi <- pnorm(s, log.p=TRUE)
  logphi <- dnorm(s, log=TRUE)
  g <- exp(logphi - logPhi)
  bad <- (!is.finite(g)) | (s < -8)
  if(any(bad)){
    t <- -s[bad]
    g[bad] <- t + 1/t
  }
  h <- g*(g + s)
  bad2 <- (!is.finite(h)) | (s < -8)
  if(any(bad2)){
    t <- -s[bad2]
    h[bad2] <- 1 + 1/(t^2)
  }
  list(g=g, h=h)
}