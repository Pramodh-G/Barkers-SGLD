# Simulate Barker's and SGLD on a chi-sq distribution to get an idea of truncated behavior.
source("utils/normal_utils.R")
library("coda")

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
N <- 1e4
x <- seq(0, 30, 0.1)
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

#--------------------- generating figure for multiple h

hs <- c(0.4, 1.2, 2, 2.8, 3.6)

chain_sgbds <- list()
ess_sgbds <- list()
oob_sgbds <- list()
chain_sglds <- list()
ess_sglds <- list()
oob_sglds <- list()

i <- 1
s_val <- 4.0
for (h in hs) {
    print("doing for h = ")
    print(h)
    chain_sgbd <- sgbd(dof, N, h, s_val)
    chain_sgbds[i] <- list(chain_sgbd)
    ess_sgbds[i] <- effectiveSize(chain_sgbd)
    oob_sgbds[i] <- sum(chain_sgbd < 0)
    print("finished sgbd ")
    chain_sgld <- sgld(dof, N, h, s_val)
    chain_sglds[i] <- list(chain_sgld)
    ess_sglds[i] <- effectiveSize(chain_sgld)
    oob_sglds[i] <- sum(chain_sgld < 0)
    print("finished sgld ")
    i <- i + 1
}

# for (i in 1:5) {
#    chain_sgld <- chain_sglds[[i]]
#    chain_sgbd <- chain_sgbds[[i]]
#    grad_sgld <- grad_log_chisq(chain_sgld)
#    grad_sgbd <- grad_log_chisq(chain_sgbd)
#    dim(chain_sgld) <- c(length(chain_sgld), 1)
#    dim(grad_sgld) <- c(length(grad_sgld), 1)
#    dim(chain_sgbd) <- c(length(chain_sgbd), 1)
#    dim(grad_sgbd) <- c(length(grad_sgbd), 1)
#    print(i)
#    print("doing for sgld")
#    ksd_sglds[i] <- ksd(chain_sgld, grad_sgld)
#    print("doing for sgbd")
#    ksd_sgbds[i] <- ksd(chain_sgbd, grad_sgbd)
# }

par(mfrow = c(1, 5))

for (i in 1:5) {
    plot(x, dchisq(x, df = dof), main = paste("desities for h:", hs[i]), col = "red", "l")
    lines(density(chain_sglds[[i]]), col = "blue")
    lines(density(chain_sgbds[[i]]), col = "green")
    legend(x = "top", c("normal density", "sgbd", "sgld"), fill = c("red", "green", "blue"))
}

png("plots/chisq-oob-h.png")
plot(x = hs, y = as.numeric(oob_sglds), "b", col="red", xlab="Step Size", ylab = "Log Kernel Stein Density")
lines(x = hs, y = as.numeric(oob_sgbds), col = "green","b")
legend(x = "topleft", c("SGBD", "SGLD"), fill = c("green", "red"))
dev.off()
