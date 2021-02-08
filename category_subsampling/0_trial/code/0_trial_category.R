library(MASS)
source("D:\\GitHub\\sub-sampling\\category_subsampling\\0_trial\\code\\functions_ss_plot.R")

#################################
#### 1 categorical variable, 2 levels:
## code 1: 1/ 0
## code 2: 1/ -1

set.seed(123)
n <- 10000
nBeta <- 7

nBeta_cont <- nBeta-2
beta_cont <- rep(0.5, nBeta_cont)
Sigma_cont <- matrix(0.5, nrow = nBeta_cont, ncol = nBeta_cont)
diag(Sigma_cont) <- 1
mu_cont <- rep(0, nBeta_cont)
X_cont <- mvrnorm(n, mu_cont, Sigma_cont)

## code 1: 1/ 0
beta0_1 <- 0.5
beta_cat1 <- 0.5
X1 <- cbind(1, X_cont, c(rep(1, n/2), rep(0, n/2)))
## code 2: 1/ -1
X2 <- cbind(1, X_cont, c(rep(1, n/2), rep(-1, n/2)))


eta <- X1 %*% c(beta0_1, beta_cont, beta_cat1)
p <- exp(eta)/(1 + exp(eta))
y <- rbinom(n, 1, prob = p)

##############################################
#### X1: 1/0 code; X2: 1/-1 code
ss_mses1 <- ss_core(y, X1)
ss_mses2 <- ss_core(y, X2)

save.image("category_ss_2levels.RData")

# subsmaple_plot(ss_mses1)
# subsmaple_plot(ss_mses2)

png('D:\\GitHub\\sub-sampling\\category_subsampling\\0_trial\\plot_2levels.png',
    width = 600,height = 600, res = 100)
compPlot(ss_mses1, ss_mses2, '1/0', '1/-1', ylim = c(0, 0.4))
dev.off()

#################################
#### 1 categorical variable, 3 levels:
## code 1: (1, 0)/ (0, 1)/ (0, 0)
## code 2: (1, 0)/ (0, 1)/ (-1, -1)
rm(list = setdiff(ls(), lsf.str()))
set.seed(234)
n <- 10000
nBeta <- 7

nBeta_cont <- nBeta-3
beta_cont <- rep(0.5, nBeta_cont)
Sigma_cont <- matrix(0.5, nrow = nBeta_cont, ncol = nBeta_cont)
diag(Sigma_cont) <- 1
mu_cont <- rep(0, nBeta_cont)
X_cont <- mvrnorm(n, mu_cont, Sigma_cont)

## code 1: 1/ 0
beta0_1 <- 0.5
beta_cat1 <- 0.5
beta_cat2 <- 0.5

## code 1: (1, 0)/ (0, 1)/ (0, 0)
X1 <- cbind(1, X_cont,
            c(rep(1, round(n/3)), rep(0, n - round(n/3))),
            c(rep(0, round(n/3)), rep(1, round(n/3)), rep(0, n - 2*round(n/3))))
## code 2: (1, 0)/ (0, 1)/ (-1, -1)
X2 <- cbind(1, X_cont,
            c(rep(1, round(n/3)), rep(0, round(n/3)), rep(-1, n - 2*round(n/3))),
            c(rep(0, round(n/3)), rep(1, round(n/3)), rep(-1, n - 2*round(n/3))))


eta <- X1 %*% c(beta0_1, beta_cont, beta_cat1, beta_cat2)
p <- exp(eta)/(1 + exp(eta))
y <- rbinom(n, 1, prob = p)

##############################################
#### X1: 1/0 code; X2: 1/-1 code
ss_mses1 <- ss_core(y, X1)
ss_mses2 <- ss_core(y, X2)

save.image("category_ss_3level.RData")

# subsmaple_plot(ss_mses1)
# subsmaple_plot(ss_mses2)

png('D:\\GitHub\\sub-sampling\\category_subsampling\\0_trial\\plot_3levels.png',
    width = 600,height = 600, res = 100)
compPlot(ss_mses1, ss_mses2, '10/01/01', '10/01/-1-1', ylim = c(0, 0.5))
dev.off()

















