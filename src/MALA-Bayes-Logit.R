# A bayesian logistic regression using the titanic dataset and MALA proposals.
# using a vector prop.sd(diagonal covariance proposal) leads to the standard MALA algorithm which does not
# allow varying levels of step sizes.

# Here, we use a preconditioning matrix (covariance matrix after a pilot run) M to allow for different step sizes.

set.seed(42)
library(mvtnorm)


titanic <- read.csv("https://dvats.github.io/assets/titanic.csv")
y <- titanic[, 1]
X <- as.matrix(titanic[, -1])

log_posterior <- function(beta)
{
    temp <- (1 - y) * X

    -sum(beta^2) / 2 - sum(log(1 + exp(-X %*% beta))) - sum(temp %*% beta)
}

grad_log_posterior <- function(beta)
{
    temp <- (1 - y) * X
    prod <- X %*% beta
    denom <- as.vector(1 / (1 + exp(prod)))
    -beta - colSums(temp) + colSums(X * denom)
}
gradient_step <- function(beta, h, M)
{
    beta + h^2 * M %*% grad_log_posterior(beta) / 2
}
mala <- function(y, X, N = 1e4, h = 0.35, M)
{
    p <- dim(X)[2]

    foo <- glm(y ~ X - 1, family = binomial("logit"))$coef
    beta <- as.matrix(foo, ncol = 1)
    beta.mat <- matrix(0, nrow = N, ncol = p)
    beta.mat[1, ] <- as.numeric(beta)

   # print(dim(grad_log_posterior(beta)))
   # print(dim(beta))
   # print(dim(rmvnorm(1, mean = numeric(p), sigma = eps)))
    accept <- 0

    for (i in 2:N)
    {
        prop <-  gradient_step(beta, h, M) + h * t(rmvnorm(1, mean = numeric(p), sigma = M))
        # prop <- beta + (eps * grad_log_posterior(beta)) / 2 + rnorm(p, mean = 0, sd = eps)
        alpha <- log_posterior(prop) - log_posterior(beta) + dmvnorm(t(beta), mean =gradient_step(prop, h, M), sigma = h^2 * M, log = TRUE) - dmvnorm(t(prop), mean = gradient_step(beta, h, M), sigma = h^2 * M, log = TRUE)


        if (log(runif(1)) < alpha)
        {
            beta <- prop
            accept <- accept + 1
        }

        beta.mat[i, ] <- beta
    }

    print(paste("acceptance: ", accept / N))
    return(beta.mat)
}

pilot <- mala(y, X, N = 1e5, h = 3.5e-3, M = diag(dim(X)[2]))
sigma_est <- cov(pilot)
chain <- mala(y, X, N = 1e5, h = 1.2, M = sigma_est)

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
