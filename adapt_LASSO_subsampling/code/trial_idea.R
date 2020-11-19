library(glmnet)
library(MASS)

#################################
#### data generation: mzNormal
set.seed(123)
n <- 10000
nBeta <- 40
beta <- rep(0, nBeta)
beta[1:10] <- 0.5
Sigma <- matrix(0.5, nrow = nBeta, ncol = nBeta)
diag(Sigma) <- 1
mu <- rep(0, nBeta)
X <- mvrnorm(n, mu, Sigma)
eta <- X %*% beta
p <- exp(eta)/(1 + exp(eta))
y <- rbinom(n, 1, prob = p)

beta.mle <- glm(y ~ X - 1,
                family = binomial(link = 'logit'))$coefficients

#################################
##### functions
adapt.lasso <- function(X, y){
  ridge1_cv <- cv.glmnet(x = X, y = y, intercept = F,
                         family = 'binomial',
                         type.measure = "mse",
                         nfold = 10, alpha = 0)
  best_ridge_coef <- as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.1se))[-1]
  
  alasso1_cv <- cv.glmnet(x = X, y = y, intercept = F,
                          family = 'binomial',
                          type.measure = "mse",
                          nfold = 10,
                          alpha = 1,
                          penalty.factor = 1 / abs(best_ridge_coef),
                          keep = TRUE)
  return(coef(alasso1_cv, s = alasso1_cv$lambda.1se)[-1])
}

#################################
r0 <- 200
r <- c(100, 200, 300, 500, 700, 1000)
S <- 1000

#### GLM not converge for small subsamples
#### adaptive LASSO
beta.adpLasso <- adapt.lasso(X, y)

## (a) full data
full.beta.boot <- matrix(NA, nrow = S, ncol = nBeta)
for(i in 1:S){
  tmp.idx <- sample(1:n, n, replace = T)
  tmp.y <- y[tmp.idx]
  tmp.X <- X[tmp.idx, ]
  full.beta.boot[i, ] <- adapt.lasso(tmp.X, tmp.y)
}

full.mse.adlasso <- mean(apply((full.beta.boot -
                          matrix(rep(beta.adpLasso, S),nrow = S, byrow = T))^2,
                       1, sum))
full.mse.true <- mean(apply((full.beta.boot -
                               matrix(rep(beta, S),nrow = S, byrow = T))^2,
                            1, sum))
full.mse.mle <- mean(apply((full.beta.boot -
                               matrix(rep(beta.mle, S),nrow = S, byrow = T))^2,
                            1, sum))

## (b) mMSE
mMSE.mse.adlasso <- rep(NA, length(r))
mMSE.mse.true <- rep(NA, length(r))
mMSE.mse.mle <- rep(NA, length(r))
for(i in 1:length(r)){
  mMSE.beta.boot <- matrix(NA, nrow = S, ncol = nBeta)
  for(j in 1:S){
    n1 <- sum(y)
    n0 <- n - n1
    pilot.ssp <- rep(1/(2*n0), n)
    pilot.ssp[y==1] <- 1/(2*n1)
    pilot.idx <- sample(1:n, r0, replace = T, prob = pilot.ssp)
    pilot.y <- y[pilot.idx]
    pilot.X <- X[pilot.idx, ]
    pilot.beta <- adapt.lasso(pilot.X, pilot.y)
    pilot.eta <- X %*% pilot.beta
    pilot.prob <- exp(pilot.eta)/(1 + exp(pilot.eta))
    pilot.prob.sub <- pilot.prob[pilot.idx]
    pilot.W <- solve(t(pilot.X)%*% (pilot.X* pilot.prob.sub*(1-pilot.prob.sub)))
    
    mMSE.ssp <- abs(y - pilot.prob)*sqrt(apply((X %*% pilot.W)^2, 1, sum))
    mMSE.ssp <- mMSE.ssp/sum(mMSE.ssp)
    mMSE.idx <- sample(1:n, r[i], replace = T, prob = mMSE.ssp)
    mMSE.y <- y[mMSE.idx]
    mMSE.X <- X[mMSE.idx, ]
    mMSE.beta <- adapt.lasso(mMSE.X, mMSE.y)
    mMSE.beta.boot[j, ] <- mMSE.beta + pilot.beta
  }
  
  mMSE.mse.adlasso[i] <- mean(apply((mMSE.beta.boot -
                               matrix(rep(beta.adpLasso, S),nrow = S, byrow = T))^2,
                            1, sum))
  mMSE.mse.true[i] <- mean(apply((mMSE.beta.boot -
                                    matrix(rep(beta, S),nrow = S, byrow = T))^2,
                                 1, sum))
  mMSE.mse.mle[i] <- mean(apply((mMSE.beta.boot -
                                    matrix(rep(beta.mle, S),nrow = S, byrow = T))^2,
                                 1, sum))
}

## (c) mVC
mVC.mse.adlasso <- rep(NA, length(r))
mVC.mse.true <- rep(NA, length(r))
mVC.mse.mle <- rep(NA, length(r))
for(i in 1:length(r)){
  mVC.beta.boot <- matrix(NA, nrow = S, ncol = nBeta)
  for(j in 1:S){
    n1 <- sum(y)
    n0 <- n - n1
    pilot.ssp <- rep(1/(2*n0), n)
    pilot.ssp[y==1] <- 1/(2*n1)
    pilot.idx <- sample(1:n, r0, replace = T, prob = pilot.ssp)
    pilot.y <- y[pilot.idx]
    pilot.X <- X[pilot.idx, ]
    pilot.beta <- adapt.lasso(pilot.X, pilot.y)
    pilot.eta <- X %*% pilot.beta
    pilot.prob <- exp(pilot.eta)/(1 + exp(pilot.eta))
    
    mVC.ssp <- abs(y - pilot.prob)*sqrt(apply(X^2, 1, sum))
    mVC.ssp <- mVC.ssp/sum(mVC.ssp)
    mVC.idx <- sample(1:n, r[i], replace = T, prob = mVC.ssp)
    mVC.y <- y[mVC.idx]
    mVC.X <- X[mVC.idx, ]
    mVC.beta <- adapt.lasso(mVC.X, mVC.y)
    mVC.beta.boot[j, ] <- mVC.beta + pilot.beta
  }
  
  mVC.mse.adlasso[i] <- mean(apply((mVC.beta.boot -
                              matrix(rep(beta.adpLasso, S),nrow = S, byrow = T))^2,
                           1, sum))
  mVC.mse.true[i] <- mean(apply((mVC.beta.boot -
                                   matrix(rep(beta, S),nrow = S, byrow = T))^2,
                                1, sum))
  mVC.mse.mle[i] <- mean(apply((mVC.beta.boot -
                                   matrix(rep(beta.mle, S),nrow = S, byrow = T))^2,
                                1, sum))
}

save.image('result.RData')

##################################################

png('D:\\GitHub\\sub-sampling\\adapt_LASSO_subsampling\\code\\MSE_adLASSO.png',
    width = 600,height = 600, res = 100)
plot(r, mMSE.mse.adlasso, type = 'b', ylim = c(0, 2),
     pch = '1', col = 1, lwd = 2,
     main = 'mzNormal_MSE_adLASSO', ylab = 'MSE_adLASSO')
lines(r, mVC.mse.adlasso, type = 'b', pch = '2',
      col = 2, lwd = 2)
abline(h = full.mse.adlasso, lty = 2, lwd = 2)
legend('topright', lwd = 2,
       col = c(1, 2, 1), pch = c('1', '2', NA),
       lty = c(1, 1, 2),
       legend = c('mMSE', 'mVc', 'full'))
dev.off()


png('D:\\GitHub\\sub-sampling\\adapt_LASSO_subsampling\\code\\MSE_MLE.png',
    width = 600,height = 600, res = 100)
plot(r, mMSE.mse.mle, type = 'b', ylim = c(0, 2),
     pch = '1', col = 1, lwd = 2,
     main = 'mzNormal_MSE_MLE', ylab = 'MSE_MLE')
lines(r, mVC.mse.mle, type = 'b', pch = '2',
      col = 2, lwd = 2)
abline(h = full.mse.mle, lty = 2, lwd = 2)
legend('topright', lwd = 2,
       col = c(1, 2, 1), pch = c('1', '2', NA),
       lty = c(1, 1, 2),
       legend = c('mMSE', 'mVc', 'full'))
dev.off()

png('D:\\GitHub\\sub-sampling\\adapt_LASSO_subsampling\\code\\MSE_true.png',
    width = 600,height = 600, res = 100)
plot(r, mMSE.mse.true, type = 'b', ylim = c(0, 2),
     pch = '1', col = 1, lwd = 2,
     main = 'mzNormal_MSE_true', ylab = 'MSE_true')
lines(r, mVC.mse.true, type = 'b', pch = '2',
      col = 2, lwd = 2)
abline(h = full.mse.true, lty = 2, lwd = 2)
legend('topright', lwd = 2,
       col = c(1, 2, 1), pch = c('1', '2', NA),
       lty = c(1, 1, 2),
       legend = c('mMSE', 'mVc', 'full'))
dev.off()














