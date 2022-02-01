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

grad_log_chisq <- function(x, dof) # pdf ∝ x^(dof/2 - 1) exp(- x / 2)
{
    if (x <= 0) {
        return(0)
    }
    - 0.5 + ((dof / 2) - 1) * (1 / x)
}

sgld <- function(dof, N, h, s_val)
{
    samples <- numeric(N)
    samples[1] <- s_val

    for (i in 2:N) {
        x <- samples[i - 1]
        samples[i] <- x + (h / 2) * grad_log_chisq(x, dof) + rnorm(1) * sqrt(h)
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
        p <- 1 / (1 + exp(-z * grad_log_chisq(x, dof)))
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

sum(samp_sgld < 0)
sum(samp_sgbd < 0)

vis_density(samp_sgbd, dof)
vis_acf(samp_sgbd)
vis_chain(samp_sgbd)

vis_density(samp_sgld, dof)
vis_acf(samp_sgld)
vis_chain(samp_sgld)
