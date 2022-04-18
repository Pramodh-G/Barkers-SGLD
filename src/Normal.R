source("utils/normal_utils.R")
library("coda")
set.seed(42)

mu <- 0
sd <- 1

N <- 1e4
x <- seq(-5 * sd + mu, 5 * sd + mu, 0.1)

# BARKER METHODS
# result <- barker_normal(N = N, h = 0.5, mu = mu, sd = sd, dist = "normal")
# chain_barker <- result$chain
# accept_barker <- result$accept
# print(accept_barker)

# result <- barker_adap_normal(N = N, mu = mu, sd = sd, dist = "normal")
# chain_barker <- result$chain
# accept_barker <- result$accept
# mu_barker <- result$mu
# sig_barker <- result$sig
# lambda_barker <- result$lambda
# print(accept_barker)
# print(mu_barker)
# print(sig_barker)
# print(lambda_barker)

chain_sgbd <- sgbd_normal(N = N, h = 0.5, mu = mu, sd = sd, dist = "normal")

# Langevin Methods
# result <- mala_normal(N = N, h = 3.4, mu = mu, sd = sd, dist = "normal")
# chain_mala <- result$chain
# accept_mala <- result$accept
# print(accept_mala)

chain_sgld <- sgld_normal(N = N, h = 2, mu = mu, sd = sd, dist = "normal")

# chain_sgld <- chain_sgbd
# acf(chain_sgld)
# plot.ts(chain_sgld)
plot(x, dnorm(x, mean = mu, sd = sd), col = "red", "l")
lines(density(chain_sgld), col = "blue")
lines(density(chain_sgbd), col = "green")

#MH Methods
# result <- mh_normal(N = N, h = 2.5, mu = mu, sd = sd, dist = "normal")
# chain_mh <- result$chain
# accept_mh <- result$accept
# print(accept_mh)

# Load Pre saved chain, mean = 0, sd = 1
chain_mh_load <- readRDS("variables/normal-mh.rds")

#--------------------- generating figure for multiple h.

hs <- c(0.2, 1, 1.6, 1.8, 2)

chain_sgbds <- list()
ess_sgbds <- list()
ksd_sgbds <- list()
chain_sglds <- list()
ess_sglds <- list()
ksd_sglds <- list()

i <- 1
for (h in hs) {
    print("doing for h = ")
    print(h)
    chain_sgbd <- sgbd_normal(N = N, h = h, mu = mu, sd = sd, dist = "normal")
    chain_sgbds[i] <- list(chain_sgbd)
    ess_sgbds[i] <- effectiveSize(chain_sgbd)
    print("finished sgbd ")
    chain_sgld <- sgld_normal(N = N, h = h, mu = mu, sd = sd, dist = "normal")
    chain_sglds[i] <- list(chain_sgld)
    ess_sglds[i] <- effectiveSize(chain_sgld)
    print("finished sgld ")
    i <- i + 1
}

for (i in 1:5) {
   chain_sgld <- chain_sglds[[i]]
   dim(chain_sgld) <- c(1, length(chain_sgld))
   chain_sgbd <- chain_sgbds[[i]]
   dim(chain_sgbd) <- c(1, length(chain_sgbd))
   ksd_sglds[i] <- ksd(chain_sgld, grad_log_norm(chain_sgld))
   ksd_sgbds[i] <- ksd(chain_sgbd, grad_log_norm(chain_sgbd))
}

# par(mfrow = c(1, 5))

# for (i in 1:5) {
#     plot(x, dnorm(x, mean = mu, sd = sd), main = paste("desities for h:", hs[i]), col = "red", "l")
#     lines(density(chain_sglds[[i]]), col = "blue")
#     lines(density(chain_sgbds[[i]]), col = "green")
#     legend(x = "top", c("normal density", "sgbd", "sgld"), fill = c("red", "green", "blue"))
# }

png("plots/normal-ksd-h.png")
plot(x = hs, y = log(as.numeric(ksd_sglds), base = 10), "b", col="red", xlab="Step Size", ylab = "Log Kernel Stein Density")
lines(x = hs, y = log(as.numeric(ksd_sgbds), base = 10), col = "green","b")
legend(x = "topleft", c("SGBD", "SGLD"), fill = c("green", "red"))
dev.off()
