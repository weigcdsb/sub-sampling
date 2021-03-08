# remove.packages('COMPoissonReg')
# .rs.restartR()
# setwd("C:\\Users\\gaw19004\\Documents\\GitHub\\sub-sampling\\COM-Poisson\\cdoe")
# install.packages("COMPoissonReg_0.7.0.tar.gz", repos = NULL, type="source")
library(COMPoissonReg)

set.seed(123)
#### generate data
n <- 10000
x <- rnorm(n, 0, 0.5)
beta <- c(0, 2)
lam <- exp(cbind(1, x) %*% beta)

g <- rnorm(n, 1, 1.5)
gam <- c(1, 0.5)
nu <- exp(cbind(1, g) %*% gam)

# plot(lam)
# plot(nu)
# abline(h = 1)
# plot(qcmp(rep(0.5, n), lam, nu))

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

# fit3: fit nu: log(nu) = g %*% gam
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

png('D:\\GitHub\\sub-sampling\\\\COM-Poisson\\cdoe\\mean.png',
    width = 600,height = 600, res = 100)
plot(y, col = 1, pch = 1, main = 'mean')
lines(fitted1, col = 2, pch = 2, type = 'p')
lines(fitted2, col = 3, pch = 3, type = 'p')
lines(fitted3, col = 4, pch = 4, type = 'p')
legend('topleft', legend = c('obs.', 'Poisson',
                             'CMP-constant nu',
                             'CMP-model nu'),
       pch = 1:4, col = 1:4)
dev.off()


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

png('D:\\GitHub\\sub-sampling\\\\COM-Poisson\\cdoe\\median.png',
    width = 600,height = 600, res = 100)
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
dev.off()


mean((y - qpois(0.5, fit1$fitted.values))^2)
mean((y - predict.cmp.quantile(0.5, fit2, cbind(1, x)))^2)
mean((y - predict.cmp.quantile(0.5, fit3, cbind(1, x)))^2)

################################################################
################################################################
#### necessary to do sub-sampling...


#################################
#### common functions
moment_sum <- function(lam, nu, nMax = 1000, gradOnly){
  
  
  if(lam >= 2 & nu <= 1){
    alph <- lam^(1/nu)
    E_y <- alph - (nu - 1)/(2*nu) -
      (nu^2 - 1)/(alph*24*(nu^2)) -
      (nu^2 - 1)/(alph^2*24*(nu^3))
    
    E_logyFac <- alph*(log(lam)/nu - 1) +
      log(lam)/(2*nu^2) + 1/(2*nu) + log(2*pi)/2 -
      (1/(24*alph))*(1 + 1/(nu^2) + log(lam)/nu - log(lam)/(nu^3)) -
      (1/(24*alph^2))*(1/(nu^3) + log(lam)/(nu^2) - log(lam)/nu^4)
    
    if(gradOnly){
      Var_y <- NA
      Var_logyFac <- NA
      Cov_y_logyFac <- NA
    }else{
      Var_y <- alph/nu + (nu^2 - 1)/(alph*24*(nu^3)) +
        (nu^2 - 1)/(alph^2*12*(nu^4))
      
      
      Var_logyFac <- alph*(log(lam))^2/(nu^3) + log(lam)/(nu^3) + 1/(2*nu^2) +
        (1/(alph*24*(nu^5)))*(-2*(nu^2) + 4*nu*log(lam) + (-1 + nu^2)*(log(lam))^2) +
        (1/((alph^2)*24*(nu^6)))*(-3*nu^2 - 2*nu*(-3 + nu^2)*log(lam) + 2*(-1 + nu^2)*(log(lam))^2)
      
      
      Cov_y_logyFac <- alph*log(lam)/(nu^2) + 1/(2*nu^2) +
        (1/(24*alph))*(2/(nu^3) + log(lam)/(nu^2) - log(lam)/(nu^4)) -
        (1/(24*alph^2))*(1/(nu^2)- 3/(nu^4) - 2*log(lam)/(nu^3) + 2*log(lam)/(nu^5))
    }
    
  }else{
    pmf <- dcmp(0:nMax, lam, nu)
    E_y <- sum((0:nMax)*pmf)
    E_logyFac <- sum((lgamma(1:(nMax + 1)))*pmf)
    
    if(gradOnly){
      Var_y <- NA
      Var_logyFac <- NA
      Cov_y_logyFac <- NA
    }else{
      Var_y <- sum((0:nMax)^2*pmf) - E_y^2
      Var_logyFac <- sum((lgamma(1:(nMax + 1))^2)*pmf) - E_logyFac^2
      Cov_y_logyFac <- sum((0:nMax)*(lgamma(1:(nMax + 1)))*pmf) -E_y*E_logyFac
    }
  }
  
  return(list(E_y = E_y, Var_y = Var_y,
              E_logyFac = E_logyFac, Var_logyFac = Var_logyFac,
              Cov_y_logyFac = Cov_y_logyFac))
}

gradInfo <- function(x, g, y, beta, gam, gradOnly = F){
  X <- cbind(1, x)
  G <- cbind(1, g)
  lam <- c(exp(X%*% beta))
  nu <- c(exp(G %*% gam))
  
  if(gradOnly){
    moments <- mapply(moment_sum, lam, nu, gradOnly = T)
    info <- NA
  }else{
    moments <- mapply(moment_sum, lam, nu, gradOnly = F)
  }
  
  grad1_each <- (y - unlist(moments['E_y', ]))* X
  grad2_each <- nu*(unlist(moments['E_logyFac', ]) - lgamma(y+1))* G
  grad_each <- cbind(grad1_each, grad2_each)
  
  grad <- apply(grad_each, 2, sum)
  
  
  if(!gradOnly){
    info1 <- t(X) %*% (unlist(moments['Var_y', ]) * X)
    info2 <- t(X) %*% (-nu*unlist(moments['Cov_y_logyFac', ])* G)
    info3 <- t(info2)
    info4 <- t(G) %*% (nu*(nu*unlist(moments['Var_logyFac', ]) -
                             unlist(moments['E_logyFac', ]) +
                             lgamma(y+1)) * G)
    
    info <- rbind(cbind(info1, info2), cbind(info3, info4))
    
  }
  
  return(list(grad_each = grad_each,
              grad = grad,
              info = info))
}


pilot <- function(x, g, y, r0){
  
  n <- length(y)
  pilot.ssp <- rep(1/n, n)
  pilot.idx <- sample(1:n, r0, replace = T, prob = pilot.ssp)
  pilot.y <- y[pilot.idx]
  pilot.x <- x[pilot.idx]
  pilot.g <- g[pilot.idx]
  pilot.ssp.star <- rep(1/r0, r0)
  
  w <- 1/pilot.ssp.star
  w <- r0*w/sum(w)
  
  fit.pilot <- glm.cmp(pilot.y ~ pilot.x, pilot.y ~ pilot.g,
                       weights = w)
  
  return(list(theta = coef(fit.pilot),
              ssp = pilot.ssp.star,
              idx = pilot.idx,
              x = pilot.x,
              g = pilot.g,
              y = pilot.y))
}

############################
S <- 500
fit.mle <- glm.cmp(y ~ x, y ~ g)
theta.mle <- coef(fit.mle)

## (a) full data
full.theta.boot <- matrix(NA, nrow = S, ncol = length(beta) + length(gam))
for(i in 1:S){
  tmp.idx <- sample(1:n, n, replace = T)
  tmp.y <- y[tmp.idx]
  tmp.x <- x[tmp.idx]
  tmp.g <- g[tmp.idx]
  
  full.theta.boot[i, ] <- coef(glm.cmp(tmp.y ~ tmp.x, tmp.y ~ tmp.g))
}

full.mse.beta <- mean(apply((full.theta.boot[, 1:length(beta)] -
                               matrix(rep(theta.mle[1:length(beta)], S),nrow = S, byrow = T))^2,
                            1, sum))
full.mse.gam <- mean(apply((full.theta.boot[, (length(beta)+1):length(theta.mle)] -
                              matrix(rep(theta.mle[(length(beta)+1):length(theta.mle)], S),nrow = S, byrow = T))^2,
                           1, sum))

full.mse <- mean(apply((full.theta.boot -
                          matrix(rep(theta.mle, S),nrow = S, byrow = T))^2,
                       1, sum))

## (b) different SSPs
r0 <- 200
# r <- c(100, 200, 300, 500, 700, 1000)
r <- c(100, 200, 300, 500, 700, 1000)
# b.1 uniform

unif.mse.beta <- rep(NA, length(r))
unif.mse.gam <- rep(NA, length(r))
unif.mse <- rep(NA, length(r))

for(i in 1:length(r)){
  unif.r.total <- r[i] + r0
  unif.theta.boot <- matrix(NA, nrow = S, ncol = length(beta) + length(gam))
  
  for(j in 1:S){
    boolFalse<-F
    while(boolFalse==F)
    {
      tryCatch({
        print(j)
        unif.idx <- sample(1:n, unif.r.total, replace = T)
        unif.y <- y[unif.idx]
        unif.x <- x[unif.idx]
        unif.g <- g[unif.idx]
        
        unif.ssp <- 1/unif.r.total
        w <- rep(1/unif.ssp, unif.r.total)
        w <- unif.r.total*w/sum(w)
        
        
        unif.theta.boot[j, ] <- coef(glm.cmp(unif.y ~ unif.x, unif.y ~ unif.g,
                                             weights = w))
        boolFalse<-T
      },error=function(e){
      },finally={})
    }
    
  }
  
  unif.mse.beta[i] <- mean(apply((unif.theta.boot[, 1:length(beta)] -
                                    matrix(rep(theta.mle[1:length(beta)], S),nrow = S, byrow = T))^2,
                                 1, sum))
  unif.mse.gam[i] <- mean(apply((unif.theta.boot[, (length(beta)+1):length(theta.mle)] -
                                   matrix(rep(theta.mle[(length(beta)+1):length(theta.mle)], S),nrow = S, byrow = T))^2,
                                1, sum))
  
  unif.mse[i] <- mean(apply((unif.theta.boot -
                               matrix(rep(theta.mle, S),nrow = S, byrow = T))^2,
                            1, sum))
}

plot(r, unif.mse.beta, type = 'l')
plot(r, unif.mse.gam, type = 'l')
plot(r, unif.mse, type = 'l')

# b.2 mMSE
mMSE.mse.beta <- rep(NA, length(r))
mMSE.mse.gam <- rep(NA, length(r))
mMSE.mse <- rep(NA, length(r))

for(i in 1:length(r)){
  
  mMSE.theta.boot <- matrix(NA, nrow = S, ncol = length(beta) + length(gam))
  
  for(j in 1:S){
    
    boolFalse<-F
    while(boolFalse==F)
    {
      tryCatch({
        print(j)
        pilot.result <- pilot(x, g, y, r0)
        
        pilot.theta <- pilot.result$theta
        gradInf <- gradInfo(x, g, y, pilot.theta[1:length(beta)],
                            pilot.theta[(length(beta) + 1):length(pilot.theta)])
        
        mMSE.ssp <- apply(abs(solve(gradInf$info) %*%
                                t(gradInf$grad_each)), 2, sum)
        mMSE.ssp <- mMSE.ssp/sum(mMSE.ssp)
        mMSE.idx <- sample(1:n, r[i], replace = T, prob = mMSE.ssp)
        
        # combined
        mMSE.y <- y[c(mMSE.idx, pilot.result$idx)]
        mMSE.x <- x[c(mMSE.idx, pilot.result$idx)]
        mMSE.g <- g[c(mMSE.idx, pilot.result$idx)]
        mMSE.ssp.star <- c(mMSE.ssp[mMSE.idx], pilot.result$ssp)
        
        w <- 1/mMSE.ssp.star
        w <- (r[i]+r0)*w/sum(w)
        
        mMSE.theta.boot[j, ] <- coef(glm.cmp(mMSE.y ~ mMSE.x, mMSE.y ~ mMSE.g,
                                             weights = w))
        boolFalse<-T
      },error=function(e){
      },finally={})
    }
    
  }
  
  mMSE.mse.beta[i] <- mean(apply((mMSE.theta.boot[, 1:length(beta)] -
                                    matrix(rep(theta.mle[1:length(beta)], S),nrow = S, byrow = T))^2,
                                 1, sum))
  mMSE.mse.gam[i] <- mean(apply((mMSE.theta.boot[, (length(beta)+1):length(theta.mle)] -
                                   matrix(rep(theta.mle[(length(beta)+1):length(theta.mle)], S),nrow = S, byrow = T))^2,
                                1, sum))
  
  mMSE.mse[i] <- mean(apply((mMSE.theta.boot -
                               matrix(rep(theta.mle, S),nrow = S, byrow = T))^2,
                            1, sum))
  
  
}

plot(r, mMSE.mse.beta, type = 'l')
plot(r, mMSE.mse.gam, type = 'l')
plot(r, mMSE.mse, type = 'l')

# b.3 mVC
mVC.mse.beta <- rep(NA, length(r))
mVC.mse.gam <- rep(NA, length(r))
mVC.mse <- rep(NA, length(r))


for(i in 1:length(r)){
  
  mVC.theta.boot <- matrix(NA, nrow = S, ncol = length(beta) + length(gam))
  
  for(j in 1:S){
    
    boolFalse<-F
    while(boolFalse==F)
    {
      tryCatch({
        print(j)
        pilot.result <- pilot(x, g, y, r0)
        
        pilot.theta <- pilot.result$theta
        grad <- gradInfo(x, g, y, pilot.theta[1:length(beta)],
                         pilot.theta[(length(beta) + 1):length(pilot.theta)], gradOnly = T)
        
        mVC.ssp <- apply(abs(grad$grad_each), 1, sum)
        mVC.ssp <- mVC.ssp/sum(mVC.ssp)
        mVC.idx <- sample(1:n, r[i], replace = T, prob = mVC.ssp)
        
        # combined
        mVC.y <- y[c(mVC.idx, pilot.result$idx)]
        mVC.x <- x[c(mVC.idx, pilot.result$idx)]
        mVC.g <- g[c(mVC.idx, pilot.result$idx)]
        mVC.ssp.star <- c(mVC.ssp[mVC.idx], pilot.result$ssp)
        
        w <- 1/mVC.ssp.star
        w <- (r[i]+r0)*w/sum(w)
        mVC.theta.boot[j, ] <- coef(glm.cmp(mVC.y ~ mVC.x, mVC.y ~ mVC.g,
                                            weights = w))
        boolFalse<-T
      },error=function(e){
      },finally={})
    }
    
  }
  
  mVC.mse.beta[i] <- mean(apply((mVC.theta.boot[, 1:length(beta)] -
                                   matrix(rep(theta.mle[1:length(beta)], S),nrow = S, byrow = T))^2,
                                1, sum))
  mVC.mse.gam[i] <- mean(apply((mVC.theta.boot[, (length(beta)+1):length(theta.mle)] -
                                  matrix(rep(theta.mle[(length(beta)+1):length(theta.mle)], S),nrow = S, byrow = T))^2,
                               1, sum))
  
  mVC.mse[i] <- mean(apply((mVC.theta.boot -
                              matrix(rep(theta.mle, S),nrow = S, byrow = T))^2,
                           1, sum))
  
}

plot(r, mVC.mse.beta, type = 'l', ylim = c(0, 1))
plot(r, mVC.mse.gam, type = 'l', ylim = c(0, 1))
plot(r, mVC.mse, type = 'l', ylim = c(0, 1))


png('D:\\GitHub\\sub-sampling\\\\COM-Poisson\\cdoe\\theta.png',
    width = 600,height = 600, res = 100)
plot(r, unif.mse, type = 'b', ylim = c(0, 0.35),
     col = 1, pch = "1", lwd = 2,
     main = 'theta', ylab = 'MSE')
lines(r, mMSE.mse, type = 'b', col = 2, pch = "2", lwd = 2)
lines(r, mVC.mse, type = 'b', col = 3, pch = "3", lwd = 2)
abline(h = full.mse, lty = 2, lwd = 2, col = 1)
legend('topright', lwd = 2,
       lty = c(rep(1, 3), 2), col = c(1:3, 1),
       pch = c(as.character(1:3), NA),
       legend = c('uniform', 'mMSE', 'mVc', 'full'))
dev.off()

png('D:\\GitHub\\sub-sampling\\\\COM-Poisson\\cdoe\\beta.png',
    width = 600,height = 600, res = 100)
plot(r, unif.mse.beta, type = 'b', ylim = c(0, 0.3),
     col = 1, pch = "1", lwd = 2,
     main = 'beta', ylab = 'MSE')
lines(r, mMSE.mse.beta, type = 'b', col = 2, pch = "2", lwd = 2)
lines(r, mVC.mse.beta, type = 'b', col = 3, pch = "3", lwd = 2)
abline(h = full.mse.beta, lty = 2, lwd = 2, col = 1)
legend('topright', lwd = 2,
       lty = c(rep(1, 3), 2), col = c(1:3, 1),
       pch = c(as.character(1:3), NA),
       legend = c('uniform', 'mMSE', 'mVc', 'full'))
dev.off()

png('D:\\GitHub\\sub-sampling\\\\COM-Poisson\\cdoe\\gamma.png',
    width = 600,height = 600, res = 100)
plot(r, unif.mse.gam, type = 'b', ylim = c(0, 0.07),
     col = 1, pch = "1", lwd = 2,
     main = 'gamma', ylab = 'MSE')
lines(r, mMSE.mse.gam, type = 'b', col = 2, pch = "2", lwd = 2)
lines(r, mVC.mse.gam, type = 'b', col = 3, pch = "3", lwd = 2)
abline(h = full.mse.gam, lty = 2, lwd = 2, col = 1)
legend('topright', lwd = 2,
       lty = c(rep(1, 3), 2), col = c(1:3, 1),
       pch = c(as.character(1:3), NA),
       legend = c('uniform', 'mMSE', 'mVc', 'full'))
dev.off()

save.image(file = 'comResults.RData')
























