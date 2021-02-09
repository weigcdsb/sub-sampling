rm(list = setdiff(ls(), lsf.str()))
set.seed(234)
n <- 2000
nBeta <- 4

nBeta_cont <- nBeta-3
beta_cont <- rep(0.5, nBeta_cont)
Sigma_cont <- matrix(0.5, nrow = nBeta_cont, ncol = nBeta_cont)
diag(Sigma_cont) <- 0.5
mu_cont <- rep(0, nBeta_cont)
X_cont <- mvrnorm(n, mu_cont, Sigma_cont)

## code 1: 1/ 0
beta0_1 <- 0.5
beta_cat1 <- 0.5
beta_cat2 <- 0.5

## code 1: (1, 0)/ (0, 1)/ (0, 0)
X1 <- cbind(1, X_cont,
            c(rep(1, round(n/4)), rep(0, n - round(n/4))),
            c(rep(0, round(n/4)), rep(1, round(n/4)), rep(0, n - 2*round(n/4))))
## code 2: (1, 0)/ (0, 1)/ (-1, -1)
X2 <- cbind(1, X_cont,
            c(rep(1, round(n/4)), rep(0, round(n/4)), rep(-1, n - 2*round(n/4))),
            c(rep(0, round(n/4)), rep(1, round(n/4)), rep(-1, n - 2*round(n/4))))


eta <- X1 %*% c(beta0_1, beta_cont, beta_cat1, beta_cat2)
p <- exp(eta)/(1 + exp(eta))
y <- rbinom(n, 1, prob = p)

ss_mses1 <- ss_core(y, X1, r = c(100, 200, 300, 500))
ss_mses2 <- ss_core(y, X2, r = c(100, 200, 300, 500))

save.image("category_ss_2levels_v2.RData")

png('D:\\GitHub\\sub-sampling\\category_subsampling\\0_trial\\plot_2levels.png',
    width = 600,height = 600, res = 100)
compPlot(ss_mses1, ss_mses2, '10/01/00', '10/01/-1-1', ylim = c(0,0.08))
dev.off()





