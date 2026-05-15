
library(actuar)

h1.fun <- function(n){
  X <- runif(n)
  U <- runif(n)
  D <- 1*(0.5*X^2 + 0.25 >= U)
  S <- rt(n, df=3)
  Y1 <- 5*S*(1+X)
  Y0 <- S*(1+X)

  res <- list(Y1=Y1, Y0=Y0, D=D, X=X)
  return(res)
}

h2.fun <- function(n){
  X <- runif(n)
  U <- runif(n)
  D <- 1*(0.5*X^2 + 0.25 >= U)
  C1 <- rfrechet(n, shape=2)
  C2 <- rfrechet(n, shape=3)
  Y1 <- C1*exp(X)
  Y0 <- C2*exp(X)

  res <- list(Y1=Y1, Y0=Y0, D=D, X=X)
  return(res)
}

h3.fun <- function(n){
  X <- runif(n)
  U <- runif(n)
  D <- 1*(0.5*X^2 + 0.25 >= U)
  Y1 <- actuar::rpareto(n, shape=1.75+X, scale=2)
  Y0 <- actuar::rpareto(n, shape=1.75+5*X, scale=1)

  res <- list(Y1=Y1, Y0=Y0, D=D, X=X)
  return(res)
}
