
library(actuar)

h1.light.fun <- function(n){
  X <- runif(n)
  U <- runif(n)
  D <- 1 * (0.5 * X^2 + 0.25 >= U)
  S <- rnorm(n)
  Y1 <- 5 * S * (1 + X)
  Y0 <- S * (1 + X)
  
  res <- list(Y1 = Y1, Y0 = Y0, D = D, X = X)
  return(res)
}


h2.light.fun <- function(n){
  X <- runif(n)
  U <- runif(n)
  D <- 1 * (0.5 * X^2 + 0.25 >= U)
  C1 <- rexp(n, rate = 1)
  C2 <- rexp(n, rate = 2)
  Y1 <- C1 * exp(X)
  Y0 <- C2 * exp(X)
  
  res <- list(Y1 = Y1, Y0 = Y0, D = D, X = X)
  return(res)
}


h3.light.fun <- function(n){
  X <- runif(n)
  U <- runif(n)
  D <- 1 * (0.5 * X^2 + 0.25 >= U)
  
  shape1 <- 2 + X
  shape0 <- 3 + 2 * X
  Y1 <- rweibull(n, shape = shape1, scale = 2)
  Y0 <- rweibull(n, shape = shape0, scale = 1)
  
  res <- list(Y1 = Y1, Y0 = Y0, D = D, X = X)
  return(res)
}

