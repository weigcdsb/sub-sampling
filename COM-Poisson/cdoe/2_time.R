
procTime <- function(n, seed = 123){
  set.seed(123)
  x <- rnorm(n, 0, 1)
  beta <- c(0, 1)
  lam <- exp(cbind(1, x) %*% beta)
  
  g <- rnorm(n, 1, 1)
  gam <- c(1, 1)
  nu <- exp(cbind(1, g) %*% gam)
  y <- rcmp(n, lam, nu)
  
  ptm <- proc.time()
  fit <- glm.cmp(formula.lambda = y ~ x, formula.nu = y ~ g)
  proc.cmp <- proc.time() - ptm
  return(proc.cmp)
}

B <- 1000
tGrid <- seq(1000, 10000, 1000)
TBoot <- matrix(NA, nrow = B, ncol = length(tGrid))
for(i in 1:B){
  print(i)
  t <- sapply(tGrid, procTime)
  TBoot[i, ] <- t[3, ]
}

























