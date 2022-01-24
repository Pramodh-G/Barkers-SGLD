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

barker <- function(y, X, N = 1e3, prop.sd = 0.35)
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
        # prop <-  gradient_step(beta, h, M) + h * t(rmvnorm(1, mean = numeric(p), sigma = M))
        grad_beta <- grad_log_posterior(beta)

        z <- rnorm(p, mean = 0, sd = prop.sd)
        prob_invert <- 1 / (1 + exp(-z * grad_log_posterior(beta)))
        inv_or_not <- (runif(p) < prob_invert)
        b <- 2 * inv_or_not - 1
        prop <- beta + z * b

        grad_prop <- grad_log_posterior(prop)

        alpha <- log_posterior(prop) - log_posterior(beta) - sum(log1p(exp( grad_prop * (prop - beta) ))) + sum(log1p(exp( grad_beta * (beta - prop) )))

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

chain <- barker(y, X, N = 1e5, prop.sd = 4.5e-3)
sigma_est <- diag(cov(chain))
chain <- barker(y, X, N = 1e5, prop.sd = 5 * sigma_est)
# nasty surprises on increasing N. Maybe needs more samples to explore space. All of a sudden drop in acceptance.


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
