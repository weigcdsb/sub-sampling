weighted.MLE <- function(X, G, y, ssp, maxIter = 1000){
  
  beta <- rep(0, ncol(X))
  gam <- rep(0, ncol(G))
  
  update <- Inf
  iter <- 0
  while((sum(update^2) > 1e-6)& (iter < maxIter)){
    
    lam <- exp(X %*% beta)
    nu <- exp(G %*% gam)
    if(iter == 0){
      lam <- rep(mean(y), length(y))
    }
    
    moments <- mapply(moment_sum, lam, nu)
    
    OP1 <- t(X) %*% (unlist(moments['Var_y', ]) * X/ ssp)
    OP2 <- t(X) %*% (-c(nu*unlist(moments['Cov_y_logyFac', ]))* G/ ssp)
    OP3 <- t(OP2)
    OP4 <- t(G) %*% (c(nu*(nu*unlist(moments['Var_logyFac', ]) -
                             unlist(moments['E_logyFac', ]) +
                             lgamma(y+1))) * G/ ssp)
    # OP4 <- t(G) %*% ((c(nu)^2)*unlist(moments['Var_logyFac', ])* G/ ssp)
    
    
    
    OP <- rbind(cbind(OP1, OP2), cbind(OP3, OP4))
    
    grad1 <- apply((y - unlist(moments['E_y', ]))* X/ ssp, 2, sum)
    grad2 <- apply(c(nu)*(unlist(moments['E_logyFac', ]) - lgamma(y+1))* G/ ssp, 2, sum)
    grad <- c(grad1, grad2)
    
    update <- solve(OP) %*% grad
    beta <- beta + update[1:ncol(X)]
    gam <- gam + update[(ncol(X) + 1): length(update)]
    
    iter <- iter + 1
    print(iter)
    print(c(beta, gam))
    print(sum(update^2))
  }
  
  if(iter < maxIter){
    return(c(beta, gam))
  }else{
    print('Not Converge')
    return()
  }
}
