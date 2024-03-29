set.seed(42)
library(mvtnorm)

# SGLD for a single dimensional normal distribution
mu  <- 3
sd <- 2

# iters
N <- 1e4
h <- 0.4
h_bim <- 1
samples <- numeric(N)
# samples_bim <- numeric(N)
samples[1] <- 4.0
# samples_bim[1] <- 4.0


for (i in 2:N) {
    samples[i] <- samples[i - 1] + h * -(samples[i - 1] - mu) / (2 * sd^2) + rnorm(1)*sqrt(h)
}

# for (i in 2:N) {
    # z <- rnorm(1, mean = sqrt(1 - h_bim), sd = sqrt(h_bim))
    # if(runif(1) < 0.5)
    # {
        # z <- -z
    # }
    # samples_bim[i] <- samples_bim[i - 1] - h_bim * samples_bim[i - 1] / 2 + z
# }

x <- seq(-5*sd + mu, 5 *sd + mu, 0.01)
plot(x, dnorm(x, mean = mu, sd = sd), type="l", col="red")
lines(density(samples))

plot(x, dnorm(x), type="l", col="red")
lines(density(samples_bim))
plot.ts(samples_bim)

# A bayesian logistic regression using the titanic dataset and MALA proposals.
# using a vector prop.sd(diagonal covariance proposal) leads to the standard MALA algorithm which does not
# allow varying levels of step sizes.

# Here, we use a preconditioning matrix (covariance matrix after a pilot run) M to allow for different step sizes.


titanic <- read.csv("https://dvats.github.io/assets/titanic.csv")
y <- titanic[, 1]
X <- as.matrix(titanic[, -1])

log_posterior <- function(beta)
{
    temp <- (1 - y) * X

    -sum(beta^2) / 2 - sum(log(1 + exp(-X %*% beta))) - sum(temp %*% beta)
}

grad_log_posterior <- function(beta, X, y)
{
    temp <- (1 - y) * X
    prod <- X %*% beta
    denom <- as.vector(1 / (1 + exp(prod)))
    -beta - colSums(temp) + colSums(X * denom)
}

sgld <- function(y, X, N = 1e4, h = 0.35, M, minibatch_size = 300)
{
    p <- dim(X)[2]
    dataset_size <- nrow(X)

    foo <- glm(y ~ X - 1, family = binomial("logit"))$coef
    beta <- as.matrix(foo, ncol = 1)
    beta.mat <- matrix(0, nrow = N, ncol = p)
    beta.mat[1, ] <- as.numeric(beta)


    for (i in 2:N)
    {
        rand_perm <- sample(dataset_size)[1:minibatch_size]
        X_mini <- X[rand_perm, ]
        y_mini <- y[rand_perm]
    beta <- beta + (dataset_size * h^2 * M %*% grad_log_posterior(beta, X_mini, y_mini)) / (2 * minibatch_size) + h * t(rmvnorm(1, mean = numeric(p), sigma = M))
        beta.mat[i, ] <- beta
    }

    return(beta.mat)
}

pilot <- sgld(y, X, N = 1e5, h = 2.9e-3, M = diag(dim(X)[2]))
sigma_est <- cov(pilot)
chain <- sgld(y, X, N = 1e5, h = 0.85, M = sigma_est)

# diagnostics.
plot.ts(chain)

par(mfrow = c(2, 3))

for (i in 1:dim(chain)[2])
{
    acf(chain[, i], main = paste("ACF for component ", i))
}

for (i in 1:dim(chain)[2])
{
    plot(density(chain[, i]), main = paste("density for component ", i))
}
