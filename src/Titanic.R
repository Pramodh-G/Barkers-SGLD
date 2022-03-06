source("utils/titanic_utils.R")
source("utils/adap_titanic_utils.R")
set.seed(42)

N <- 1e4
p <- 6
prior_sig <- 5

# Barker methods
# result <- barker_titanic(y, X, N = N, h = 4e-3, dist = "normal")
# chain_barker <- result$chain
# accept_barker <- result$accept
# print(accept_barker)
# chain_sgbd <- sgbd_titanic(y, X, N = N, h = 1e-2, dist = "normal")

# Langevin methods
# chain_sgld <- sgld_titanic(y, X, N = N, h = 2.5e-3, dist = "normal")

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
# chain_mh_load <- readRDS("variables/titanic-mh.rds")

# chain <- chain_adap_mh

plot.ts(chain)

par(mfrow = c(2, 3))

for (i in 1:p)
{
    acf(chain[, i], main = paste("ACF for component ", i))
    # plot(result$mus[, i], main = paste("mean for component ", i))
    # plot(result$sig[, i], main = paste("mean for component ", i))
}

for (i in 1:dim(chain)[2])
{
    plot(density(chain[, i]), main = paste("density for component ", i))
}
