source("utils/common_utils.R")

log_norm <- function(x, mu = 0, sd = 1)
{
    ret <- ((-(x - mu)^2) / (2 * sd^2))
    return(ret)
}

grad_log_norm <- function(x, mu = 0, sd = 1)
{
    ret <- ((-(x - mu)) / sd^2)
    return(ret)
}

######### BARKER ###################
barker_adap_normal <- function(N = 1e3, mu = 0, sd = 1, dist = "normal")
{
    samples <- numeric(N)
    samples[1] <- 4.0
    accept <- 0

    lambda <- (2.4^2)
    mu_adap <- 0
    sigma_adap <- 1

    for (i in 2:N) {
        x <- samples[i - 1]
        z <- samp_z(n = 1, h = lambda * sigma_adap, dist = dist)
        p <- 1 / (1 + exp(-z * grad_log_norm(x, mu = mu, sd = sd)))
        b <- -1
        if(runif(1) < p)
        {
            b <- 1
        }
        prop <- x + b * z

        log_alpha <- log_norm(prop, mu, sd) - log_norm(x, mu, sd) - log1p(exp( (prop - x) * grad_log_norm(prop, mu = mu, sd = sd) )) + log1p(exp( (x - prop) * grad_log_norm(x, mu =mu, sd =sd) ))
        if(log(runif(1)) < log_alpha)
        {
            samples[i] <- prop
            accept <- accept + 1
        }
        else
        {
            samples[i] <- x
        }
        gamma <- (i^(-0.6))
        lambda <- exp(log(lambda) + gamma * (min(1, exp(log_alpha)) - 0.44))
        sigma_adap <- sigma_adap + gamma * ((samples[i] - mu_adap)^2 - sigma_adap)
        mu_adap <- mu_adap + gamma * (samples[i] - mu_adap)
    }
    ret <- list(chain = samples, accept = accept / N, mu = mu_adap, sig = sigma_adap, lambda = lambda)
    return(ret)
}
# dist can be "normal" or "bim_normal"
barker_normal <- function(N = 1e3, h = 0.2, mu = 0, sd = 1, dist = "normal")
{
    samples <- numeric(N)
    samples[1] <- 4.0
    accept <- 0

    for (i in 2:N) {
        x <- samples[i - 1]
        z <- samp_z(n = 1, h = h, dist = dist)
        p <- 1 / (1 + exp(-z * grad_log_norm(x, mu = mu, sd = sd)))
        b <- -1
        if(runif(1) < p)
        {
            b <- 1
        }
        prop <- x + b * z

        log_alpha <- log_norm(prop, mu, sd) - log_norm(x, mu, sd) - log1p(exp( (prop - x) * grad_log_norm(prop, mu = mu, sd = sd) )) + log1p(exp( (x - prop) * grad_log_norm(x, mu =mu, sd =sd) ))
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
    ret <- list(chain = samples, accept = accept / N)
    return(ret)
}

# dist can be "normal", "bim_normal"
sgbd_normal <- function(N = 1e3, h = 0.2, mu = 0, sd = 1, dist = "normal")
{
    samples <- numeric(N)
    samples[1] <- 4.0

    for (i in 2:N) {
        x <- samples[i - 1]
        z <- samp_z(n = 1, h = h, dist = "normal")
        p <- 1 / (1 + exp(-z * grad_log_norm(x, mu = mu, sd = sd)))
        b <- -1
        if(runif(1) < p)
        {
            b <- 1
        }
        samples[i] <- x + b * z
    }
    return(samples)
}

######### LANGEVIN ##########
mala_normal <- function(N = 1e3, h = 0.2, mu = 0, sd = 1, dist = "normal")
{
    samples <- numeric(N)
    samples[1] <- 4.0
    accept <- 0
    for (i in 2:N) {
       z <- samp_z(n = 1, h = h, dist = dist)
       x <- samples[i - 1]
       prop <- x + h^2 * grad_log_norm(x, mu = mu, sd = sd) / 2  + z
       grad_step_prop <- prop + h^2 * grad_log_norm(prop, mu = mu, sd = sd) / 2
       grad_step_cur  <- x + h^2 * grad_log_norm(x, mu = mu, sd = sd) / 2

       log_alpha <- log_norm(prop, mu = mu, sd = sd) - log_norm(x, mu = mu, sd = sd)
       log_alpha <- log_alpha + log_proposal_dist(x, mean = grad_step_prop, sd = h)
       log_alpha <- log_alpha - log_proposal_dist(prop, mean = grad_step_cur, sd = h)

       if (log(runif(1)) < log_alpha) {
           accept <- accept + 1
          samples[i] <- prop
       }
       else {
          samples[i] <- samples[i - 1]
       }
    }
    ret <- list(chain = samples, accept = accept / N)
    return(ret)
}

sgld_normal <- function(N = 1e3, h = 0.2, mu = 0, sd = 1, dist = "normal")
{
    samples <- numeric(N)
    samples[1] <- 4.0

    for (i in 2:N) {
        z <- samp_z(n = 1, h = h, dist = dist)
        samples[i] <- samples[i - 1] + h^2 * grad_log_norm(samples[i - 1], mu = mu, sd = sd) / 2 + z
    }
    return(samples)
}

########## METROPOLIS-HASTINGS
mh_normal <- function(N = 1e3, h = 0.2, mu = 0, sd = 1, dist = "normal")
{
    samples <- numeric(N)
    samples[1] <- 4.0
    accept <- 0

    for (i in 2:N) {
        z <- samp_z(n = 1, h = h, dist = dist)
        samples[i] <- samples[i - 1] + z

        log_alpha <- log_norm(samples[i], mu = mu, sd = sd) - log_norm(samples[i - 1], mu = mu, sd = sd)

        if(log(runif(1)) < log_alpha)
        {
            accept <- accept + 1
        }
        else {
           samples[i] <- samples[i - 1]
        }
    }
    ret <- list(chain = samples, accept = accept / N)
    return(ret)
}