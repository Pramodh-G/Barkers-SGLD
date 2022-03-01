source("utils/common_utils.R")
library("mvtnorm")

# initialize it in the original file :)
prior_sig <- NULL

log_posterior <- function(beta, X, y)
{
    temp <- (1 - y) * X

    -sum(beta^2) / (2 * prior_sig^2) - sum(log(1 + exp(-X %*% beta))) - sum(temp %*% beta)
}

grad_log_posterior <- function(beta, X, y)
{
    temp <- (1 - y) * X
    prod <- X %*% beta
    denom <- as.vector(1 / (1 + exp(prod)))
    (-beta / (prior_sig^2)) - colSums(temp) + colSums(X * denom)
}

gradient_step <- function(beta, X, y, h)
{
    beta + h^2 * grad_log_posterior(beta, X, y) / 2
}

######### Barker
barker_aryt <- function(y, X, N = 1e3, h = 0.6, dist = "normal")
{
    p <- dim(X)[2]

    foo <- glm(y ~ X - 1, family = binomial("logit"))$coef
    beta <- as.matrix(foo, ncol = 1)
    col_na <- which(is.na(beta))
    beta[col_na] <- colMeans(X)[col_na]
    beta.mat <- matrix(0, nrow = N, ncol = p)
    beta.mat[1, ] <- as.numeric(beta)

    accept <- 0

    for (i in 2:N)
    {
        grad_beta <- grad_log_posterior(beta, X, y)

        z <- samp_z(n = p, h = h, dist = dist)
        prob_invert <- 1 / (1 + exp(-z * grad_beta))
        inv_or_not <- (runif(p) < prob_invert)
        b <- 2 * inv_or_not - 1
        prop <- beta + z * b

        grad_prop <- grad_log_posterior(prop, X, y)

        log_alpha <- log_posterior(prop, X, y) - log_posterior(beta, X, y) - sum(log1p(exp( grad_prop * (prop - beta) ))) + sum(log1p(exp( grad_beta * (beta - prop) )))

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

sgbd_aryt <- function(y, X, N = 1e3, h = 0.35, minibatch_size = 452, dist = "normal")
{
    p <- dim(X)[2]
    dataset_size <- nrow(X)

    foo <- glm(y ~ X - 1, family = binomial("logit"))$coef
    beta <- as.matrix(foo, ncol = 1)
    col_na <- which(is.na(beta))
    beta[col_na] <- colMeans(X)[col_na]
    beta.mat <- matrix(0, nrow = N, ncol = p)
    beta.mat[1, ] <- as.numeric(beta)


    for (i in 2:N)
    {
        rand_perm <- sample(dataset_size)[1:minibatch_size]
        X_mini <- X[rand_perm, ]
        y_mini <- y[rand_perm]
        z <- samp_z(n = p, h = h, dist = dist)
        prob_invert <- 1 / (1 + exp(-z * grad_log_posterior(beta, X_mini, y_mini)))
        inv_or_not <- (runif(p) < prob_invert)
        b <- 2 * inv_or_not - 1
        beta.mat[i, ] <- beta + z * b
        beta <- beta.mat[i, ]
    }

    return(beta.mat)
}

##### LANGEVIN

sgld_aryt <- function(y, X, N = 1e3, h = 0.35, minibatch_size = 452, dist = "normal")
{
    p <- dim(X)[2]
    dataset_size <- nrow(X)
    foo <- glm(y ~ X - 1, family = binomial("logit"))$coef
    beta <- as.matrix(foo, ncol = 1)
    col_na <- which(is.na(beta))
    beta[col_na] <- colMeans(X)[col_na]
    beta.mat <- matrix(0, nrow = N, ncol = p)
    beta.mat[1, ] <- as.numeric(beta)


    for (i in 2:N)
    {
        rand_perm <- sample(dataset_size)[1:minibatch_size]
        X_mini <- X[rand_perm, ]
        y_mini <- y[rand_perm]
        z <- samp_z(n = p, h = h, dist = "normal")
        beta <- beta + h^2 * grad_log_posterior(beta, X_mini, y_mini) / 2 + z
        beta.mat[i, ] <- beta
    }

    return(beta.mat)
}

mala_aryt <- function(y, X, N = 1e4, h = 0.35)
{
    p <- dim(X)[2]

    foo <- glm(y ~ X - 1, family = binomial("logit"))$coef
    beta <- as.matrix(foo, ncol = 1)
    col_na <- which(is.na(beta))
    beta[col_na] <- colMeans(X)[col_na]
    beta.mat <- matrix(0, nrow = N, ncol = p)
    beta.mat[1, ] <- as.numeric(beta)

    accept <- 0

    for (i in 2:N)
    {
        M <- diag(p)
        prop <-  gradient_step(beta, X, y, h) + h * t(rmvnorm(1, mean = numeric(p), sigma = h^2 * M))
        alpha <- log_posterior(prop, X, y) - log_posterior(beta, X, y) + dmvnorm(t(beta), mean =gradient_step(prop, X, y, h), sigma = h^2 * M, log = TRUE) - dmvnorm(t(prop), mean = gradient_step(beta, X, y, h), sigma = h^2 * M, log = TRUE)


        if (log(runif(1)) < alpha)
        {
            beta <- prop
            accept <- accept + 1
        }

        beta.mat[i, ] <- beta
    }

    ret <- list(chain = beta.mat, accept = accept / N)
    return(ret)
}

###### METROPOLIS- HASTINGS
mh_aryt <- function(y, X, N = 1e4, h = 0.35, dist = "normal")
{
    p <- dim(X)[2]

    model_glm <- glm(y ~ X - 1, family = binomial("logit"))
    foo <- model_glm$coef
    beta <- as.matrix(foo, ncol = 1)
    col_na <- which(is.na(beta))
    beta[col_na] <- colMeans(X)[col_na]
    beta.mat <- matrix(0, nrow = N, ncol = p)
    beta.mat[1, ] <- as.numeric(beta)
    accept <- 0

    for (i in 2:N)
    {
        z <- samp_z(n = p, h = h, dist = dist)
        prop <- beta + z
        alpha <- log_posterior(prop, X, y) - log_posterior(beta, X, y)
        if (log(runif(1)) < alpha)
        {
            beta <- prop
            accept <- accept + 1
        }

        beta.mat[i, ] <- beta
    }

    ret <- list(chain = beta.mat, accept = accept / N)
    return(ret)
}









# Code procured from authors, used to reproduce results as accurately as possible!
aryt.uci <- read.csv(url("https://archive.ics.uci.edu/ml/machine-learning-databases/arrhythmia/arrhythmia.data"),
                     header=FALSE)

### DATA PREPROCESSING ########
y<-as.numeric(aryt.uci[,280]==1)
X<-aryt.uci[,-280]
## Remove columns with missing values (denoted as "?") or identical terms
cols.to.remove<-c()
for(j in 1:dim(aryt.uci)[2]){
  if(any(aryt.uci[,j]=="?")){
    cols.to.remove<-c(cols.to.remove,j)
  }else{
    if(sd(aryt.uci[,j])==0){
      cols.to.remove<-c(cols.to.remove,j)
    }
  }
}
X<-aryt.uci[,-c(cols.to.remove,280)]

# import and select the correct columns (columns are imported because the sample.int function is unstable
# across R versions)
load(file = "data/cols.RData")
# import cols.RData before this line to get access to the 'cols' variable
X <- X[, cols]

n<-dim(X)[1]
p<-dim(X)[2]

# # rescale if necessary
X <- scale(X)
# prior_variance <- 25 # choose prior variance