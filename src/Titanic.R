source("utils/titanic_utils.R")
source("utils/adap_titanic_utils.R")
source("utils/common_utils.R")
set.seed(42)

N <- 1e4
p <- 6
prior_sig <- 5

# Barker methods
result <- barker_titanic(y, X, N = 1e4, h = 4e-3, dist = "normal")
chain_barker <- result$chain
accept_barker <- result$accept
print(accept_barker)
pre_sd <- sqrt(diag(cov(chain_barker)))

chain_sgbd <- sgbd_titanic(y, X, N = N, h = 3e-1 * pre_sd, dist = "normal")
chain_sgbd_grad <- prep_gradients(chain_sgbd)
ksd(chain_sgbd, chain_sgbd_grad) # AROUND 59.37793

# Langevin methods
chain_sgld <- sgld_titanic(y, X, N = N, h = 4e-1 * pre_sd, dist = "normal")
chain_sgld_grad <- prep_gradients(chain_sgld)
ksd(chain_sgld, chain_sgld_grad) # AROUND 28.6686

# result <- mala_titanic(y, X, N = N, h = 3.85e-3)
# chain_mala <- result$chain
# accept_mala <- result$accept
# print(accept_mala)

# MH Methods
# result <- mh_titanic(y, X, N = N, h = 4e-3, dist = "normal")
# chain_mh <- result$chain
# accept_mh <- result$accept
# print(accept_mh)

# result <- mh_adap_titanic(y, X, N = N, dist = "normal")
# plot(result$lambdas)
# chain_adap_mh <- result$chain
# print(result$accept)
# print(result$sig)
# print(result$mu)
# print(result$lambda)

# Load presaved MH
chain_mh_load <- readRDS("variables/titanic-mh.rds")

chain <- chain_mh_load

# plot.ts(chain)

par(mfrow = c(2, 3))

# for (i in 1:p)
# {
#     acf(chain[, i], main = paste("ACF for component ", i))
#     # plot(result$mus[, i], main = paste("mean for component ", i))
#     # plot(result$sig[, i], main = paste("mean for component ", i))
# }

for (i in 1:dim(chain)[2])
{
    plot(density(chain[, i]), main = paste("density for component ", i), col="black")
    lines(density(chain_sgld[, i]), col = "red")
    lines(density(chain_sgbd[, i]), col = "blue")
}
