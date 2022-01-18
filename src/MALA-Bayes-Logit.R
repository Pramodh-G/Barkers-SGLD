# A bayesian logistic regression using the titanic dataset and MALA proposals.

set.seed(42)
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
    denom <- as.vector( 1 / (1 + exp(prod)) )
    -beta - colSums(temp) + colSums(X * denom)
}

mala <- function(y, X, N = 1e4, prop.sd = 0.35)
{
    p <- dim(X)[2]

    foo <- glm(y ~ X - 1, family = binomial("logit"))$coef
    beta <- as.matrix(foo, ncol = 1)
    # print(as.numeric(beta))
    beta.mat <- matrix(0, nrow = N, ncol = p)
    # print(beta.mat)
    beta.mat[1, ] <- as.numeric(beta)
    accept <- 0

    for (i in 2:N)
    {
        prop <- beta + prop.sd^2 * grad_log_posterior(beta) + rnorm(p, mean = 0, sd = prop.sd)

        alpha <- log_posterior(prop) - log_posterior(beta)

        if (log(runif(1)) < alpha)
        {
            beta <- prop
            accept <- accept + 1
        }

        beta.mat[i, ] <- beta
    }

    print(paste("acceptance: ", accept / N))
    print(accept)
    return(beta.mat)
}

chain <- mala(y, X, N = 1e3, prop.sd = 0.009)

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