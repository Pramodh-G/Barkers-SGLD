source("utils/normal_utils.R")
set.seed(42)

mu <- 5
sd <- 3

N <- 1e6
x <- seq(-5 * sd + mu, 5 * sd + mu, 0.1)

# BARKER METHODS
# result <- barker_normal(N = N, h = 0.5, mu = mu, sd = sd, dist = "normal")
# chain_barker <- result$chain
# accept_barker <- result$accept
# print(accept_barker)

result <- barker_adap_normal(N = N, mu = mu, sd = sd, dist = "normal")
chain_barker <- result$chain
accept_barker <- result$accept
mu_barker <- result$mu
sig_barker <- result$sig
lambda_barker <- result$lambda
print(accept_barker)
print(mu_barker)
print(sig_barker)
print(lambda_barker)


# chain_sgbd <- sgbd_normal(N = N, h = 0.5, mu = mu, sd = sd, dist = "bim")

# Langevin Methods
# result <- mala_normal(N = N, h = 3.4, mu = mu, sd = sd, dist = "normal")
# chain_mala <- result$chain
# accept_mala <- result$accept
# print(accept_mala)

# chain_sgld <- sgld_normal(N = N, h = 2, mu = mu, sd = sd, dist = "laplace")

chain_sgld <- chain_barker
acf(chain_sgld)
plot.ts(chain_sgld)
plot(density(chain_sgld), col = "blue")
lines(x, dnorm(x, mean = mu, sd = sd), col = "red")

#MH Methods
# result <- mh_normal(N = N, h = 2.5, mu = mu, sd = sd, dist = "normal")
# chain_mh <- result$chain
# accept_mh <- result$accept
# print(accept_mh)

# Load Pre saved chain
chain_mh_load <- readRDS("variables/normal-mh.rds")


acf(chain_mh_load)
# plot.ts(chain_mh)
plot(density(chain_mh_load), col = "blue")
lines(x, dnorm(x, mean = mu, sd = sd), col = "red")
lines(density(chain_sgld), col = "green")
