set.seed(42)

# generate points from a N(0, 1) distribution using SGLD.
mu  <- 0
sigma <- 1

# iters
N <- 1e4
h <- 2
samples <- numeric(N)
samples[1] <- 4.0
accept <- 0

for (i in 2:N) {
    x <- samples[i - 1]
    z <- rnorm(1, mean = 0, sd = h)
    p <- 1 / (1 + exp(z * x))
    # print(p)
    b <- -1
    if(runif(1) < p)
    {
        b <- 1
    }
    prop <- x + b * z

    log_alpha <- (-prop^2 / 2) - (-samples[i - 1]^2) / 2 - log1p(exp( (prop - x) * (-prop) )) + log1p(exp( (x - prop) * (-x) ))
    if(log(runif(1)) < log_alpha)
    {
        samples[i] <- prop
        accept <- accept + 1
    }
    else
    {
       samples[i] <- x
    }
}

print(accept / N)
plot.ts(samples)
acf(samples)
x <- seq(-5, 5, 0.01)
plot(x, dnorm(x), type="l", col="blue")
lines(density(samples), col = "red")
