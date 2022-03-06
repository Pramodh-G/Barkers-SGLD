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

gradient_step <- function(beta, h)
{
    beta + h^2 * grad_log_posterior(beta) / 2
}

gamma <- function(iteration)
{
    return(((iteration)^(-0.6)))
}

upd_lambda <- function(lambda, curr_gamma, alpha, opt_alpha)
{
    log_ans <- log(lambda) + curr_gamma * (alpha - opt_alpha)
    return(exp(log_ans))
}

upd_mu <- function(mu, curr_gamma, x_next)
{
    new_mu <- mu + curr_gamma * (x_next - mu)
    return(new_mu)
}

upd_sig <- function(sig, curr_gamma, x_next, mu)
{
    vec <- x_next - mu
    # sig <- diag(sig)
    new_sigma <- sig + curr_gamma * (vec^2 - sig)
    if(!is.numeric(new_sigma)) {print("add another diag here.")}
    return(new_sigma)
}

#### Adaptive rwm? ###
mh_adap_titanic <- function(y, X, N = 1e3, dist = "normal")
{
    p <- dim(X)[2]
    opt_alpha <- 0.234

    # foo <- glm(y ~ X - 1, family = binomial("logit"))$coef
    # beta <- as.numeric(foo, ncol = 1)
    accept <- 0

    beta.mat <- matrix(0, nrow = N, ncol = p)
    mus <- matrix(0, nrow = N, ncol = p)
    sigs <- matrix(0, nrow = N, ncol = p)
    lambdas <- numeric(N)
    alphas <- numeric(N)

    beta <- numeric(p)
    mu <- rep(0, p)
    sig <- rep(1, p)
    lambda <- (2.4)^2 / (p ^ (1 / 3))
    lambdas[1] <- lambda
    beta.mat[1, ] <- beta


    for (i in 2:N)
    {
        if(!is.numeric(lambda * sig)) {throw("h is not numeric boiii")}
        z <- rnorm(n = p, sd = sqrt(lambda * sig))
        # prop <- t(rmvnorm(n = 1, mean = beta, sigma = diag(lambda * sig)))
        prop <- beta + z

        log_alpha <- log_posterior(prop) - log_posterior(beta)
        alpha <- min(1, exp(log_alpha))

        if(runif(1) < alpha)
        {
            beta <- prop
            accept <- accept + 1
        }

        curr_gamma <- gamma(i)
        lambda <- upd_lambda(lambda, curr_gamma, alpha, opt_alpha)
        mu <- upd_mu(mu, curr_gamma, beta)
        sig <- upd_sig(sig, curr_gamma, beta, mu)

        beta.mat[i, ] <- beta
        lambdas[i] <- lambda
        alphas[i] <- alpha
        sigs[i, ] <- sig
        mus[i, ] <- mu

        # print(sprintf("i: %d, l_a: %f, a: %f, l: %f", i, log_alpha, alpha, lambda))
        # print(mu)
        # print(sig)
    }

    ret <- list(chain = beta.mat, accept = accept / N, sig = sig, mus = mus, lambdas = lambdas, alphas = alphas)
    return(ret)
}