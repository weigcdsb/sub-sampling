library(COMPoissonReg)

data("freight")
mean(freight$broken)
var(freight$broken)
# over-dispersed

#### standard Poisson

fit.poi <- glm(broken ~ transfers,
               data = freight,
               family = poisson, na.action = na.exclude)
fit.poi

#### CMP
# for fair comparison, nu = constant
fit.cmp <- glm.cmp(broken ~ transfers, data = freight)
fit.cmp
coef(fit.cmp)
nu(fit.cmp)

# sdev(fit.cmp)
# vcov(fit.cmp)
# LRT: POI vs. CMP
# equitest(fit.cmp)
# deviance(fit.cmp)


#### comparison
plot(freight$broken, type = 'b', col = 1, pch = 1)
lines(fit.poi$fitted.values, type = 'b', col = 2, pch = 2)
lines(predict(fit.cmp, newdata=freight), type = 'b', col = 3, pch = 3)


new.data <- data.frame(transfers=(0:10))
y.hat <- predict(fit.cmp, newdata=new.data)
plot(0:10, y.hat, type="l",
     xlab="number of transfers", ylab="predicted number broken")

predict(fit.cmp, newdata=freight)

fit.cmp.lam <- function(fit, X){
  return(exp(as.matrix(X) %*% coef(fit)[1:ncol(fit$X)]))
}

lam <- fit.cmp.lam(fit.cmp)
nu <- nu(fit.cmp)
predict.cmp.quantile <- function(q, fit){
  
  lam <- fit.cmp.lam(fit, cbind(1, freight$transfers))
  nu <- nu(fit)
  out <- rep(NA, nrow(cbind(1, freight$transfers)))
  for(i in 1:nrow(cbind(1, freight$transfers))){
    out[i] <- qcmp(q, lam[i], nu[i])
  }
  return(out)
}

plot(predict(fit.cmp, newdata=freight), ylim = c(0, 30))
lines(predict.cmp.quantile(0.025, fit.cmp), lty = 2)
lines(predict.cmp.quantile(0.5, fit.cmp))
lines(predict.cmp.quantile(0.975, fit.cmp), lty = 2)
lines(fit.poi$fitted.values, col = 2)










