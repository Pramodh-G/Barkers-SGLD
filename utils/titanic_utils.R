source("utils/common_utils.R")
library("mvtnorm")

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

######### Barker
barker_titanic <- function(y, X, N = 1e3, h = 0.35, dist = "normal")
{
    p <- dim(X)[2]

    foo <- glm(y ~ X - 1, family = binomial("logit"))$coef
    beta <- as.matrix(foo, ncol = 1)
    beta.mat <- matrix(0, nrow = N, ncol = p)
    beta.mat[1, ] <- as.numeric(beta)

    accept <- 0

    for (i in 2:N)
    {
        grad_beta <- grad_log_posterior(beta)

        z <- samp_z(n = p, h = h, dist = dist)
        prob_invert <- 1 / (1 + exp(-z * grad_log_posterior(beta)))
        inv_or_not <- (runif(p) < prob_invert)
        b <- 2 * inv_or_not - 1
        prop <- beta + z * b

        grad_prop <- grad_log_posterior(prop)

        log_alpha <- log_posterior(prop) - log_posterior(beta) - sum(log1p(exp( grad_prop * (prop - beta) ))) + sum(log1p(exp( grad_beta * (beta - prop) )))

        if (log(runif(1)) < log_alpha)
        {
            beta <- prop
            accept <- accept + 1
        }

        beta.mat[i, ] <- beta
    }
    ret <- list(chain = beta.mat, accept = accept / N)
    return(ret)
}

sgbd_titanic <- function(y, X, N = 1e3, h = 0.35, dist = "normal")
{
    p <- dim(X)[2]

    foo <- glm(y ~ X - 1, family = binomial("logit"))$coef
    beta <- as.matrix(foo, ncol = 1)
    beta.mat <- matrix(0, nrow = N, ncol = p)
    beta.mat[1, ] <- as.numeric(beta)


    for (i in 2:N)
    {
        z <- samp_z(n = p, h = h, dist = dist)
        prob_invert <- 1 / (1 + exp(-z * grad_log_posterior(beta)))
        inv_or_not <- (runif(p) < prob_invert)
        b <- 2 * inv_or_not - 1
        beta.mat[i, ] <- beta + z * b
        beta <- beta.mat[i, i]
    }

    return(beta.mat)
}

##### LANGEVIN

sgld <- function(y, X, N = 1e3, h = 0.35, dist = "normal")
{
    p <- dim(X)[2]
    foo <- glm(y ~ X - 1, family = binomial("logit"))$coef
    beta <- as.matrix(foo, ncol = 1)
    beta.mat <- matrix(0, nrow = N, ncol = p)
    beta.mat[1, ] <- as.numeric(beta)


    for (i in 2:N)
    {
        z <- samp_z(n = p, h = h, dist = "normal")
        beta <- beta + h^2 * grad_log_posterior(beta) / 2 + z
        beta.mat[i, ] <- beta
    }

    return(beta.mat)
}

#TODO: IMPLEMENT MALA ON THIS!!!!
###### METROPOLIS- HASTINGS
mh <- function(y, X, N = 1e4, prop.sd = 0.35, dist = "normal")
{
    p <- dim(X)[2]

    foo <- glm(y ~ X - 1, family = binomial("logit"))$coef
    beta <- as.matrix(foo, ncol = 1)
    beta.mat <- matrix(0, nrow = N, ncol = p)
    beta.mat[1, ] <- as.numeric(beta)
    accept <- 0

    for (i in 2:N)
    {
        prop <- beta + samp_z(n = p, h = h, dist = dist)

        alpha <- log_posterior(prop) - log_posterior(beta)

        if (log(runif(1)) < alpha)
        {
            beta <- prop
            accept <- accept + 1
        }

        beta.mat[i, ] <- beta
    }

    ret <- list(chain = beta.mat, accept <- accept/N)
}