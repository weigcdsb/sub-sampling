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
#### different sub-sampling procedure
## uniform, A-optimization & L-optimization

ss_core <- function(y, X, S = 1000, r0 = 200,
                    r = c(100, 200, 300, 500, 700, 1000)){
  n <- dim(X)[1]
  nBeta <- dim(X)[2]
  
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
  
  full.mse <- mean(apply((X %*% t(full.beta.boot -
                                    matrix(rep(beta.mle, S),
                                           nrow = S, byrow = T)))^2,
                         2, sum))/(n - nBeta)
  
  
  ## (b) different SSPs
  
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
    
    unif.mse[i] <- mean(apply((X %*% t(unif.beta.boot -
                                         matrix(rep(beta.mle, S),
                                                nrow = S, byrow = T)))^2,
                              2, sum))/(n - nBeta)
    
  }
  
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
    
    mMSE.mse[i] <- mean(apply((X %*% t(mMSE.beta.boot -
                                         matrix(rep(beta.mle, S),
                                                nrow = S, byrow = T)))^2,
                              2, sum))/(n - nBeta)
    
  }
  
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
    
    mVC.mse[i] <- mean(apply((X %*% t(mVC.beta.boot -
                                        matrix(rep(beta.mle, S),
                                               nrow = S, byrow = T)))^2,
                             2, sum))/(n - nBeta)
    
  }
  
  return(list(r = r, 
              unif.mse = unif.mse,
              mMSE.mse = mMSE.mse,
              mVC.mse = mVC.mse,
              full.mse = full.mse))
}


#################################
## plot 1
subsmaple_plot <- function(ss_mses){
  r <- ss_mses$r
  unif.mse <- ss_mses$unif.mse
  mMSE.mse <- ss_mses$mMSE.mse
  mVC.mse <- ss_mses$mVC.mse
  full.mse <- ss_mses$full.mse
  
  plot(r, unif.mse, type = 'b', ylim = c(0, 1),
       col = 1, pch = "1", lwd = 2,
       main = 'mzNormal', ylab = 'MSE')
  lines(r, mMSE.mse, type = 'b', col = 2, pch = "2", lwd = 2)
  lines(r, mVC.mse, type = 'b', col = 3, pch = "3", lwd = 2)
  abline(h = full.mse, lty = 2, lwd = 2, col = 1)
  legend('topright', lwd = 2,
         lty = c(rep(1, 3), 2), col = c(1:3, 1),
         pch = c(as.character(1:3), NA),
         legend = c('uniform', 'mMSE', 'mVc', 'full'))
}

## compare plot
compPlot <- function(ss_mses1, ss_mses2,
                     code1, code2, ylim){
  
  plot(ss_mses1$r, ss_mses1$unif.mse, type = 'b', ylim = ylim,
       col = 1, pch = "1", lwd = 2, lty = 1,
       xlab = 'r', ylab = 'MSE', main = paste(code1, 'vs.', code2))
  lines(ss_mses2$r, ss_mses2$unif.mse,
        type = 'b', col = 1, pch = "1", lwd = 2, lty = 2)
  
  lines(ss_mses1$r, ss_mses1$mMSE.mse,
        type = 'b', col = 2, pch = "2", lwd = 2, lty = 1)
  lines(ss_mses2$r, ss_mses2$mMSE.mse,
        type = 'b', col = 2, pch = "2", lwd = 2, lty = 2)
  
  lines(ss_mses1$r, ss_mses1$mVC.mse,
        type = 'b', col = 3, pch = "3", lwd = 2, lty = 1)
  lines(ss_mses2$r, ss_mses2$mVC.mse,
        type = 'b', col = 3, pch = "3", lwd = 2, lty = 2)
  
  abline(h = ss_mses1$full.mse, lty = 1, lwd = 2, col = 4)
  abline(h = ss_mses2$full.mse, lty = 2, lwd = 2, col = 4)
  
  legend('topright', lwd = 2,
         lty = rep(1:2, 4), col = rep(1:4, rep(2, 4)),
         pch = c(as.character(rep(1:3, rep(2, 3))), NA, NA),
         legend = c(paste('uniform:', code1),
                    paste('uniform:', code2),
                    paste('mMSE:', code1),
                    paste('mMSE:', code2),
                    paste('mVc:', code1),
                    paste('mVc:', code2),
                    paste('full:', code1),
                    paste('full:', code2)))
  
}



















