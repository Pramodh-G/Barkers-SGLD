source("common_utils.R")

log_norm <- function(x, mu = 0, sd = 1)
{
    ret <- (-(x - mu)^2) / (2 * sd^2)
    return(ret)
}

grad_log_norm <- function(x, mu = 0, sd = 1)
{
    ret <- (-(x - mu)) / sd^2
    return(ret)
}

######### BARKER ###################
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

        # log_alpha <- (-(prop - mu)^2 / (2 * sigma^2)) - (-(samples[i - 1] - mu)^2) / (2 * sigma^2) - log1p(exp( (prop - x) * (-(prop - mu)/(sigma^2)) )) + log1p(exp( (x - prop) * (-(x - mu)/(sigma^2)) ))
        log_alpha <- log_norm(prop, mu, sd) - log_norm(x, mu, sd) - log1p(exp( (prop - x) * grad_log_norm(prop) )) + log1p(exp( (x - prop) * grad_log_norm(x) ))
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
mala <- function(N = 1e3, h = 0.2, mu = 0, sd = 1, dist = "normal")
{
    samples <- numeric(N)
    samples[1] <- 4.0
    accept <- 0
    for (i in 2:N) {
       z <- samp_z(n = 1, h = h, dist = dist)
       samples[i] <- samples[i - 1] + h^2 * grad_log_norm(samples[i - 1], mu = mu, sd = sd) / 2  + z

       log_alpha <- log_norm(samples[i], mu = mu, sd = sd) - log_norm(samples[i - 1], mu = mu, sd = sd)
       grad_step_prop <- samples[i] + h * grad_log_norm(samples[i], mu = mu, sd = sd) / 2
       grad_step_cur  <- samples[i - 1] + h * grad_log_norm(samples[i - 1], mu = mu, sd = sd) / 2

       log_alpha <- log_alpha + log_proposal_dist(samples[i - 1], mean = grad_step_prop, h = h, dist = dist)
       log_alpha <- log_alpha - log_proposal_dist(samples[i], mean = grad_step_cur, h = h, dist = dist)

       if (log(runif(1)) < log_alpha) {
           accept <- accept + 1
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
mh <- function(N = 1e3, h = 0.2, mu = 0, sd = 1, dist = "normal")
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