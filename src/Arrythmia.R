source("utils/arrhythmia_utils.R")

# Note that we already define X, and y in the corresponding utils section.
set.seed(42)

N <- 1e4
prior_sig <- 5

# Barker methods
#  pilot_result <- barker_aryt(y, X, N = N, h = 1.1e-1, dist = "normal")
# accept_barker <- pilot_result$accept
# print(accept_barker)
# pre_sd <- sqrt(diag(cov(pilot_result$chain)))
# result <- barker_aryt(y, X, N = N, h = 2e-1 * pre_sd, dist = "normal")
# chain_barker <- result$chain
# accept_barker <- result$accept
# print(accept_barker)
# chain_sgbd <- sgbd_aryt(y, X, N = N, h = 5e-1 * pre_sd, dist = "normal")


# Langevin methods
# pilot <- sgld_aryt(y, X, N = N, h = 0.5, dist = "normal")
# pre_sd <- sqrt(diag(cov(pilot)))
# chain_sgld <- sgld_aryt(y, X, N = 1e4, h = 0.1 * pre_sd, dist = "normal")

# result <- mala_aryt(y, X, N = 1e4, h = 0.1 * pre_sd)
# chain_mala <- result$chain
# accept_mala <- result$accept
# print(accept_mala)

# MH Methods
# pilot_result <- mh_aryt(y, X, N = N, h = 7e-2, dist = "normal")
# accept_mh <- pilot_result$accept
# print(accept_mh)
# pre_sd <- sqrt(diag(cov(pilot_result$chain)))
# result <- mh_aryt(y, X, N = N, h = 1.5e-1 * pre_sd, dist = "normal")
# chain_mh <- result$chain
# accept_mh <- result$accept
# print(accept_mh)

# Load presaved MH
chain_mh_load <- readRDS("variables/aryt-mh.rds")

chain <- chain_barker

nplot <- 10
plot.ts(chain[, 1:nplot])

par(mfrow = c(2, nplot/2))

for (i in 1:nplot)
{
    acf(chain[, i], main = paste("ACF for component ", i))
}

for (i in 1:nplot)
{
    plot(density(chain[, i]), main = paste("density for component ", i))
    lines(density(chain_mh_load[, i]), main = paste("density for component ", i, col = "blue"))
}
