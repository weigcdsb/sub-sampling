library(lars)
library(MASS)

set.seed(123)
n <- 10000
nBeta <- 50
beta <- rep(0, nBeta)
beta[1:10] <- 5
Sigma <- matrix(0.5, nrow = nBeta, ncol = nBeta)
diag(Sigma) <- 1
mu <- rep(0, nBeta)
X <- mvrnorm(n, mu, Sigma)
y <- X %*% beta + rnorm(n, 0, 15)

beta.mle <- lm(y ~ X - 1)$coefficients
# plot(beta.mle)

##############
frac <- seq(0, 1, length.out = 10)
r <- c(100, 200, 500, 1000, length(y))
S <- 1000

############# CV
# cv.index <- vector(mode = 'list', length = length(r))
# cv.cv <- vector(mode = 'list', length = length(r))
# for(i in 1:length(r)){
# 
#   tmp.idx <- sample(1:n, r[i], replace = T)
#   tmp.y <- y[tmp.idx]
#   tmp.X <- X[tmp.idx, ]
#   tmp.cv <- cv.lars(tmp.X, tmp.y,
#                     intercept = F, type = 'lasso', plot.it = F)
#   cv.index[[i]] <- tmp.cv$index
#   cv.cv[[i]] <- tmp.cv$cv
# }
# plot(cv.index[[1]], cv.cv[[1]], type = 'b', pch = '1')
# for(i in 2:length(r)){
#   lines(cv.index[[i]], cv.cv[[i]], col = i,
#         pch = as.character(i), type = 'b')
# }

##################

MSE.true <- matrix(NA, nrow = length(r), ncol = length(frac))
for(i in 1:length(r)){
  sse.true.boot <- matrix(NA, nrow = S, ncol = length(frac))
  for(j in 1:S){
    
    tmp.idx <- sample(1:n, r[i], replace = T)
    tmp.y <- y[tmp.idx]
    tmp.X <- X[tmp.idx, ]
    tmp.object <- lars(tmp.X, tmp.y,
                       intercept = F, type = 'lasso')
    tmp.coef <- coef(tmp.object, s = frac, mode = 'fraction')
    sse.true.boot[j, ] <- apply((tmp.coef -
                                   matrix(rep(beta, nrow(tmp.coef)),
                                          nrow = nrow(tmp.coef),
                                          byrow = T))^2, 1, sum)
    
  }
  MSE.true[i, ] <- apply(sse.true.boot, 2, mean)
}

########################
save.image('D:\\GitHub\\sub-sampling\\LASSO_subsampling\\2_lars_subsample\\bootResult.RData')

png('D:\\GitHub\\sub-sampling\\LASSO_subsampling\\2_lars_subsample\\MSE_true.png',
    width = 600,height = 600, res = 100)
markerPts <- seq(1, length(frac), 2)
plot(frac, MSE.true[1, ], type = 'l', lwd = 2,
     col = 1, ylim = c(0, max(MSE.true)),
     xlab = '|beta|/ max(|beta|)', ylab = 'MSE')
points(frac[markerPts], MSE.true[1, markerPts], pch = '1', col = 1)
for(i in 2:length(r)){
  
  if(i < length(r)){
    lines(frac, MSE.true[i, ], col = i, lwd = 2)
    points(frac[markerPts], MSE.true[i, markerPts],
           pch = as.character(i), col = i)
  }else{
    lines(frac, MSE.true[i, ], col = 1, lty = 2, lwd = 2)
  }
  
}
legend('topright', lwd = 2,
       col = c(1:(length(r)-1), 1),
       lty = c(rep(1, length(r)-1), 2),
       pch = c(as.character(1:(length(r)-1)), NA),
       legend = c(r[1:(length(r)-1)], 'full'))

dev.off()



