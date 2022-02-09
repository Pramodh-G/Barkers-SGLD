source("utils/normal_utils.R")
set.seed(42)

mu <- 0
sd <- 1

N <- 1e7
x <- seq(-5 * sd + mu, 5 * sd + mu, 0.1)

# BARKER METHODS
# result <- barker_normal(N = N, h = 0.5, mu = mu, sd = sd, dist = "normal")
# chain_barker <- result$chain
# accept_barker <- result$accept
# print(accept_barker)

# chain_sgbd <- sgbd_normal(N = N, h = 0.5, mu = mu, sd = sd, dist = "bim")

# Langevin Methods
# result <- mala_normal(N = N, h = 3.4, mu = mu, sd = sd, dist = "normal")
# chain_mala <- result$chain
# accept_mala <- result$accept
# print(accept_mala)

# chain_sgld <- sgld_normal(N = N, h = 1.4, mu = mu, sd = sd, dist = "normal")

# acf(chain_mala)
# plot.ts(chain_mala)
# plot(density(chain_mala), col = "blue")
# lines(x, dnorm(x, mean = mu, sd = sd), col = "red")

#MH Methods
result <- mh_normal(N = N, h = 2.5, mu = mu, sd = sd, dist = "normal")
chain_mh <- result$chain
accept_mh <- result$accept
print(accept_mh)

acf(chain_mh)
# plot.ts(chain_mh)
plot(density(chain_mh), col = "blue")
lines(x, dnorm(x, mean = mu, sd = sd), col = "red")
