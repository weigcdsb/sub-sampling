library(MASS)

#################################
#### data generation: mzNormal
set.seed(123)
n <- 10000
nBeta <- 7
beta <- rep(0.5, nBeta)
Sigma <- matrix(0.5, nrow = nBeta, ncol = nBeta)
diag(Sigma) <- 1
mu <- rep(0, nBeta)
X <- mvrnorm(n, mu, Sigma)
eta <- X %*% beta
p <- exp(eta)/(1 + exp(eta))
y <- rbinom(n, 1, prob = p)


#################################
#### common functions
weighted.MLE <- function(X, y, ssp, maxIter = 1000){
  
  beta <- rep(0, ncol(X))
  update <- Inf
  iter <- 0
  while((sum(update^2) > 1e-6)& (iter < maxIter)){
    eta <- X %*% beta
    prob <- c(exp(eta)/(1 + exp(eta)))
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

pilot <- function(X, y, r0){
  n1 <- sum(y)
  n0 <- n - n1
  pilot.ssp <- rep(1/(2*n0), n)
  pilot.ssp[y==1] <- 1/(2*n1)
  pilot.idx <- sample(1:n, r0, replace = T, prob = pilot.ssp)
  pilot.y <- y[pilot.idx]
  pilot.X <- X[pilot.idx, ]
  pilot.ssp.star <- pilot.ssp[pilot.idx]
  return(list(beta = weighted.MLE(pilot.X, pilot.y, pilot.ssp.star),
              ssp = pilot.ssp.star,
              idx = pilot.idx,
              X = pilot.X, y = pilot.y))
}

#################################
#### fit the data
S <- 1000
beta.mle <- glm(y ~ X - 1,
                family = binomial(link = 'logit'))$coefficients
## (a) full data
full.beta.boot <- matrix(NA, nrow = S, ncol = nBeta)
for(i in 1:S){
  tmp.idx <- sample(1:n, n, replace = T)
  tmp.y <- y[tmp.idx]
  tmp.X <- X[tmp.idx, ]
  full.beta.boot[i, ] <- glm(tmp.y ~ tmp.X - 1,
                             famil = binomial(link = 'logit'))$coefficients
}

full.mse <- mean(apply((full.beta.boot -
                          matrix(rep(beta.mle, S),nrow = S, byrow = T))^2,
                       1, sum))

## (b) different SSPs
r0 <- 200
r <- c(100, 200, 300, 500, 700, 1000)

# b.1 uniform

unif.mse <- rep(NA, length(r))
for(i in 1:length(r)){
  unif.r.total <- r[i] + r0
  unif.beta.boot <- matrix(NA, nrow = S, ncol = nBeta)
  
  for(j in 1:S){
    unif.idx <- sample(1:n, unif.r.total, replace = T)
    unif.y <- y[unif.idx]
    unif.X <- X[unif.idx, ]
    unif.ssp <- 1/n
    unif.beta.boot[j, ] <- t(weighted.MLE(unif.X, unif.y, unif.ssp))
  }
  
  unif.mse[i] <- mean(apply((unif.beta.boot -
                               matrix(rep(beta.mle, S),nrow = S, byrow = T))^2,
                            1, sum))
  
}

# plot(r, unif.mse, type = 'l', ylim = c(0, 1))

# b.2 mMSE
mMSE.mse <- rep(NA, length(r))
for(i in 1:length(r)){
  
  mMSE.beta.boot <- matrix(NA, nrow = S, ncol = nBeta)
  for(j in 1:S){
    pilot.result <- pilot(X, y, r0)
    
    pilot.beta <- pilot.result$beta
    pilot.eta <- X %*% pilot.beta
    pilot.prob <- exp(pilot.eta)/(1 + exp(pilot.eta))
    pilot.prob.sub <- pilot.prob[pilot.result$idx]
    pilot.W <- solve(t(pilot.result$X)%*% (pilot.result$X* 
                                             pilot.prob.sub*(1-pilot.prob.sub)/
                                             pilot.result$ssp))
    mMSE.ssp <- abs(y - pilot.prob)*sqrt(apply((X %*% pilot.W)^2, 1, sum))
    mMSE.ssp <- mMSE.ssp/sum(mMSE.ssp)
    mMSE.idx <- sample(1:n, r[i], replace = T, prob = mMSE.ssp)
    
    # combined
    mMSE.y <- y[c(mMSE.idx, pilot.result$idx)]
    mMSE.X <- X[c(mMSE.idx, pilot.result$idx), ]
    mMSE.ssp.star <- c(mMSE.ssp[mMSE.idx], pilot.result$ssp)
    mMSE.beta.boot[j, ] <- t(weighted.MLE(mMSE.X, mMSE.y, mMSE.ssp.star))
  }
  
  mMSE.mse[i] <- mean(apply((mMSE.beta.boot -
                               matrix(rep(beta.mle, S),nrow = S, byrow = T))^2,
                            1, sum))
  
}

# plot(r, unif.mse, type = 'l', ylim = c(0, 1))
# lines(r, mMSE.mse, col = 'red')

# b.3 mVC
mVC.mse <- rep(NA, length(r))
for(i in 1:length(r)){
  
  mVC.beta.boot <- matrix(NA, nrow = S, ncol = nBeta)
  for(j in 1:S){
    pilot.result <- pilot(X, y, r0)
    
    pilot.beta <- pilot.result$beta
    pilot.eta <- X %*% pilot.beta
    pilot.prob <- exp(pilot.eta)/(1 + exp(pilot.eta))
    
    mVC.ssp <- abs(y - pilot.prob)*sqrt(apply(X^2, 1, sum))
    mVC.ssp <- mVC.ssp/sum(mVC.ssp)
    mVC.idx <- sample(1:n, r[i], replace = T, prob = mVC.ssp)
    
    # combined
    mVC.y <- y[c(mVC.idx, pilot.result$idx)]
    mVC.X <- X[c(mVC.idx, pilot.result$idx), ]
    mVC.ssp.star <- c(mVC.ssp[mVC.idx], pilot.result$ssp)
    mVC.beta.boot[j, ] <- t(weighted.MLE(mVC.X, mVC.y, mVC.ssp.star))
  }
  
  mVC.mse[i] <- mean(apply((mVC.beta.boot -
                               matrix(rep(beta.mle, S),nrow = S, byrow = T))^2,
                            1, sum))
  
}

# plot(r, unif.mse, type = 'l', ylim = c(0, 1))
# lines(r, mVC.mse, col = 'red')

# b.4 LCC
LCC.mse <- rep(NA, length(r))
for(i in 1:length(r)){
  LCC.beta.boot <- matrix(NA, nrow = S, ncol = nBeta)
  for(j in 1:S){
    
    n1 <- sum(y)
    n0 <- n - n1
    pilot.ssp <- rep(1/(2*n0), n)
    pilot.ssp[y==1] <- 1/(2*n1)
    pilot.idx <- sample(1:n, r0, replace = T, prob = pilot.ssp)
    pilot.y <- y[pilot.idx]
    pilot.X <- X[pilot.idx, ]
    pilot.beta <- glm(pilot.y ~ pilot.X - 1,
                      family = binomial(link = 'logit'))$coefficients
    pilot.eta <- X %*% pilot.beta
    pilot.prob <- exp(pilot.eta)/(1 + exp(pilot.eta))
    
    LCC.ssp <- abs(y - pilot.prob)
    
    # fair comparison: subsampling with replacement
    LCC.idx <- sample(1:n, r[i], replace = T, prob = LCC.ssp)
    
    LCC.y <- y[LCC.idx]
    LCC.X <- X[LCC.idx, ]
    LCC.beta <- glm(LCC.y ~ LCC.X - 1,
                    family = binomial(link = 'logit'))$coefficients
    LCC.beta.boot[j, ] <- LCC.beta + pilot.beta
  }
  
  LCC.mse[i] <- mean(apply((LCC.beta.boot -
                              matrix(rep(beta.mle, S),nrow = S, byrow = T))^2,
                           1, sum))
  
}

# plot(r, unif.mse, type = 'l', ylim = c(0, 1))
# lines(r, LCC.mse, col = 'red')

save.image(file = "mzNormal_rep_result.RData")

#################################
#### Let's plot!
png('D:\\GitHub\\sub-sampling\\OSMAC_reproduce\\plot.png',
    width = 600,height = 600, res = 100)
plot(r, unif.mse, type = 'b', ylim = c(0, 1),
     col = 1, pch = "1", lwd = 2,
     main = 'mzNormal', ylab = 'MSE')
lines(r, mMSE.mse, type = 'b', col = 2, pch = "2", lwd = 2)
lines(r, mVC.mse, type = 'b', col = 3, pch = "3", lwd = 2)
lines(r, LCC.mse, type = 'b', col = 4, pch = "4", lwd = 2)
abline(h = full.mse, lty = 2, lwd = 2, col = 1)
legend('topright', lwd = 2,
       lty = c(rep(1, 4), 2), col = c(1:4, 1),
       pch = c(as.character(1:4), NA),
       legend = c('uniform', 'mMSE', 'mVc', 'LCC', 'full'))
dev.off()
















