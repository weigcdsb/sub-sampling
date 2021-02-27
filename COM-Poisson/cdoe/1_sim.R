library(COMPoissonReg)

set.seed(123)
#### generate data
n <- 10000
x <- rnorm(n, 0, 0.5)
beta <- c(0, 2)
lam <- exp(cbind(1, x) %*% beta)

g <- rnorm(n, -1, 1)
gamma <- c(1, 0.5)
nu <- exp(cbind(1, g) %*% gamma)

plot(lam)
plot(nu)
abline(h = 1)
plot(qcmp(rep(0.5, n), lam, nu))

y <- rcmp(n, lam, nu)

## check

# fit1: Poisson regression
ptm <- proc.time()
fit1 <- glm(y ~ x, family = poisson)
proc.poi <- proc.time() - ptm

# fit2: assume constant nu
ptm <- proc.time()
fit2 <- glm.cmp(formula.lambda = y ~ x)
proc.cmp.nuOff <- proc.time() - ptm

# fit3: fit nu: log(nu) = g %*% gamma
ptm <- proc.time()
fit3 <- glm.cmp(formula.lambda = y ~ x, formula.nu = y ~ g)
proc.cmp.nuOn <- proc.time() - ptm

## processing time
rbind(proc.poi, proc.cmp.nuOff, proc.cmp.nuOn)[, 1:3]

## coefficient
coef(fit1)
coef(fit2)
coef(fit3)

## nu
# plot(nu, col = 'red', pch = 1)
# abline(h = 1, lty = 2)
# abline(h = nu(fit2)[1])
# lines(nu(fit3), type = 'p', pch = 2)

## fitted values (mean)
fitted1 <- predict(fit1, newdata = data.frame(x = x, g = g),
                   type = 'response')
fitted2 <- predict(fit2, newdata = data.frame(x = x, g = g))
fitted3 <- predict(fit3, newdata = data.frame(x = x, g = g))


plot(y, col = 1, pch = 1, main = 'mean')
lines(fitted1, col = 2, pch = 2, type = 'p')
lines(fitted2, col = 3, pch = 3, type = 'p')
lines(fitted3, col = 4, pch = 4, type = 'p')
legend('topleft', legend = c('obs.', 'Poisson',
                             'CMP-constant nu',
                             'CMP-model nu'),
       pch = 1:4, col = 1:4)

mean((y - fitted1)^2)
mean((y - fitted2)^2)
mean((y - fitted3)^2)

## compare the median
predict.cmp.quantile <- function(q, fit, X){
  lam <- exp(X %*% coef(fit)[1:ncol(fit$X)])
  nu <- nu(fit)
  out <- rep(NA, nrow(X))
  for(i in 1:nrow(X)){
    out[i] <- qcmp(q, lam[i], nu[i])
  }
  return(out)
}

plot(y, col = 1, pch = 1, main = 'median')
lines(qpois(rep(0.5, n), fit1$fitted.values),
      type = 'p', col = 2, pch = 2)
lines(predict.cmp.quantile(0.5, fit2, cbind(1, x)),
      type = 'p', col = 3, pch = 3)
lines(predict.cmp.quantile(0.5, fit3, cbind(1, x)),
      type = 'p', col = 4, pch = 4)
legend('topleft', legend = c('obs.', 'Poisson',
                             'CMP-constant nu',
                             'CMP-model nu'),
       pch = 1:4, col = 1:4)

mean((y - qpois(0.5, fit1$fitted.values))^2)
mean((y - predict.cmp.quantile(0.5, fit2, cbind(1, x)))^2)
mean((y - predict.cmp.quantile(0.5, fit3, cbind(1, x)))^2)

################################################################
################################################################
#### necessary to do sub-sampling...


#################################
#### common functions
moment_sum <- function(lam, nu, nMax = 10000, tol = 1e-6){
  term <- rep(NA, 6)
  cum <- rep(0, 6)
  for(i in 1:nMax){
    
    term <- c(exp((i-1)*log(lam) - nu*lgamma(i)),
              exp(log(i-1) + (i-1)*log(lam) - nu*lgamma(i)),
              exp( 2*log(i-1) + (i-1)*log(lam) - nu*lgamma(i)),
              exp((i-1)*log(lam) - nu*lgamma(i))*lgamma(i),
              exp((i-1)*log(lam) - nu*lgamma(i))*(lgamma(i))^2,
              exp(log(i-1) + (i-1)*log(lam) - nu*lgamma(i))*lgamma(i))
    
    !is.na(max(pmax(term)/ pmin(term))) & 
    
    if(i > 5 & !is.na(max(pmax(term)/ pmin(term)))){
      if(max(pmax(term)/ pmin(term)) < tol){
        break
      }
    }
    cum <- cum + term
  }
  
  E_y <- cum[2]/cum[1]
  Var_y <- cum[3]/cum[1] - E_y^2
  E_logyFac <- cum[4]/cum[1]
  Var_logyFac <- cum[5]/cum[1] - E_logyFac^2
  Cov_y_logyFac <- cum[6]/cum[1] - E_y*E_logyFac
  return(list(E_y = E_y, Var_y = Var_y,
              E_logyFac = E_logyFac, Var_logyFac = Var_logyFac,
              Cov_y_logyFac = Cov_y_logyFac))
}


weighted.MLE <- function(X, G, y, ssp, maxIter = 1000){
  
  beta <- rep(0, ncol(X))
  gamma <- rep(0, ncol(G))
  
  update <- Inf
  iter <- 0
  while((sum(update^2) > 1e-6)& (iter < maxIter)){
    
    lam <- exp(X %*% beta)
    nu <- exp(G %*% gamma)
    moments <- mapply(moment_sum, lam, nu)
    
    
    
    
    
    OP <- t(X) %*% (prob * (1 - prob)* X/ ssp)
    update <- solve(OP) %*% apply((y - prob)* X/ ssp, 2, sum)
    beta <- beta + update
    iter <- iter + 1
  }
  
  
  if(iter < maxIter){
    return(beta)
  }else{
    print('Not Converge')
    return()
  }
}































