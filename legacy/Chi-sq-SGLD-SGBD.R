# Simulate Barker's and SGLD on a chi-sq distribution to get an idea of truncated behavior.

set.seed(42)

vis_acf <- function(chain)
{
    acf(chain)
}

vis_density <- function(chain, dof )
{
    x <- seq(0, 30, 0.1)
    plot(x, dchisq(x, df = dof), type = "l")
    lines(density(chain), col = "red")
}

vis_chain <- function(chain)
{
    plot.ts(chain)
}

grad_log_chisq <- function(x, dof, mode = "mala") # pdf ∝ x^(dof/2 - 1) exp(- x / 2)
{
    if (x <= 0) {
        if(mode == "mala")
        {
            return(0)
        }
        else {
           return(1e7)
        }
    }
    - 0.5 + ((dof / 2) - 1) * (1 / x)
}

sgld <- function(dof, N, h, s_val)
{
    samples <- numeric(N)
    samples[1] <- s_val

    for (i in 2:N) {
        x <- samples[i - 1]
        samples[i] <- x + (h / 2) * grad_log_chisq(x, dof, mode="mala") + rnorm(1) * sqrt(h)
    }

    return(samples)
}

#SGBD is simply removing acceptance step here, there is only one gradient
sgbd <- function(dof, N, h, s_val)
{
    samples <- numeric(N)
    samples[1] <- s_val

    for (i in 2:N) {
        x <- samples[i - 1]
        z <- rnorm(1, mean = 0, sd = h)
        p <- 1 / (1 + exp(-z * grad_log_chisq(x, dof, mode = "barker")))
        b <- -1
        if(runif(1) < p)
        {
            b <- 1
        }
        samples[i] <- x + b * z
    }

    return(samples)
}

sgbd_bim <- function(dof, N, h, s_val)
{
    samples <- numeric(N)
    samples[1] <- s_val

    for (i in 2:N) {
        x <- samples[i - 1]
        # generate samples from normals at √(1 - σ²), with variance σ²
        z <- rnorm(1, mean = (sqrt(1 - h^2)), sd = h)
        if(runif(1) < 0.5)
        {
            z <- -z
        }
        p <- 1 / (1 + exp(-z * grad_log_chisq(x, dof, mode="barkers")))
        b <- -1
        if(runif(1) < p)
        {
            b <- 1
        }
        samples[i] <- x + b * z
    }

    return(samples)
}

dof <- 10
N <- 1e5
samp_sgld <- sgld(dof, N, 0.6, 4.0)
samp_sgbd <- sgbd(dof, N, 0.9, 4.0)
samp_sgbd_bim <- sgbd_bim(dof, N, 0.9, 4.0)

sum(samp_sgld < 0)
sum(samp_sgbd < 0)
sum(samp_sgbd_bim < 0)

vis_density(samp_sgbd, dof)
vis_acf(samp_sgbd)
vis_chain(samp_sgbd)

vis_density(samp_sgld, dof)
vis_acf(samp_sgld)
vis_chain(samp_sgld)

vis_density(samp_sgbd_bim, dof)
vis_acf(samp_sgbd_bim)
vis_chain(samp_sgbd_bim)
