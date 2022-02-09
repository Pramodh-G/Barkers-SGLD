set.seed(42)
library(mvtnorm)

# generate points from a N(0, 1) distribution using barker's proposal.
mu  <- 3
sd <- 2

# iters
N <- 1e5
h <- 0.5
samples <- numeric(N)
samples[1] <- 4.0

for (i in 2:N) {
    x <- samples[i - 1]
    z <- rnorm(1, mean = 0, sd = h)
    p <- 1 / (1 + exp(-z * -(x - mu)/ (sd^2)))
    b <- -1
    if(runif(1) < p)
    {
        b <- 1
    }
    samples[i] <- x + b * z
}

plot.ts(samples)
acf(samples)
x <- seq(-5 + mu, 5 + mu, 0.01)
plot(x, dnorm(x, mean = mu, sd = sd), type ="l", col="blue")
lines(density(samples), col = "red")


# A bayesian logistic regression using the titanic dataset and barker's proposals.
# using a vector prop.sd(diagonal covariance proposal) leads to the standard barker's algorithm which does not
# allow varying levels of step sizes.
titanic <- read.csv("https://dvats.github.io/assets/titanic.csv")
y <- titanic[, 1]
X <- as.matrix(titanic[, -1])

grad_log_posterior <- function(beta)
{
    temp <- (1 - y) * X
    prod <- X %*% beta
    denom <- as.vector(1 / (1 + exp(prod)))
    -beta - colSums(temp) + colSums(X * denom)
}

barker <- function(y, X, N = 1e3, prop.sd = 0.35)
{
    p <- dim(X)[2]

    foo <- glm(y ~ X - 1, family = binomial("logit"))$coef
    beta <- as.matrix(foo, ncol = 1)
    beta.mat <- matrix(0, nrow = N, ncol = p)
    beta.mat[1, ] <- as.numeric(beta)


    for (i in 2:N)
    {
        z <- rnorm(p, mean = 0, sd = prop.sd)
        prob_invert <- 1 / (1 + exp(-z * grad_log_posterior(beta)))
        inv_or_not <- (runif(p) < prob_invert)
        b <- 2 * inv_or_not - 1
        beta.mat[i, ] <- beta + z * b
        beta <- beta.mat[i, i]
    }

    return(beta.mat)
}

barker_bim <- function(y, X, N = 1e3, prop.sd = 0.35)
{
    p <- dim(X)[2]

    foo <- glm(y ~ X - 1, family = binomial("logit"))$coef
    beta <- as.matrix(foo, ncol = 1)
    beta.mat <- matrix(0, nrow = N, ncol = p)
    beta.mat[1, ] <- as.numeric(beta)


    for (i in 2:N)
    {
        z <- rnorm(p, mean = (sqrt(1 - prop.sd^2)), sd = prop.sd)
        for(k in 1:p)
        {
            if(runif(1) < 0.5)
            {
                z[k] <- -z[k]
            }
        }

        prob_invert <- 1 / (1 + exp(-z * grad_log_posterior(beta)))
        inv_or_not <- (runif(p) < prob_invert)
        b <- 2 * inv_or_not - 1
        beta.mat[i, ] <- beta + z * b
    }

    return(beta.mat)
}

chain <- barker(y, X, N = 1e3, prop.sd = 4.5e-3)
sigma_est <- diag(cov(chain))
chain <- barker(y, X, N = 1e5, prop.sd = 4 * sigma_est)
# nasty surprises on increasing N. Maybe needs more samples to explore space. All of a sudden drop in acceptance.

chain_bim <- barker_bim(y, X, N = 1e3, prop.sd = 4.5e-3)
sigma_est_bim <- diag(cov(chain))
chain_bim <- barker(y, X, N = 1e5, prop.sd = 4 * sigma_est)


# diagnostics.
plot.ts(chain_bim)

par(mfrow = c(2, 3))

for (i in 1:dim(chain)[2])
{
    acf(chain_bim[, i], main = paste("ACF for component ", i))
}

for (i in 1:dim(chain)[2])
{
    plot(density(chain_bim[, i]), main = paste("density for component ", i))
}
