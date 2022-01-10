set.seed(42)

# generate points from a N(0, 1) distribution using SGLD.
mu  <- 0
sigma <- 1

# iters
N <- 1e4
h <- 1
samples <- numeric(N)
samples[1] <- 4.0


for (i in 2:N) {
    samples[i] <- samples[i - 1] - h * samples[i - 1] / 2 + rnorm(1)*sqrt(h)
}

hist(samples, breaks = 30)
