source("utils/common_utils.R")
library("mvtnorm")

# implement adaptive MCMC for the arrythmia dataset.

# initialize it in the original file :)
prior_sig <- NULL

#-------------------HELPER-FUNCTIONS------------------

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

gamma <- function(iteration)
{
    return((iteration)^-0.6)
}

upd_lambda <- function(curr_lambda, iteration, alpha, alpha_opt)
{
    # note alpha opt is optimal acceptance rate, 0.574 for gradient based methods, and 0.234 for MH.
    # alpha is current acceptance rate.
    log_lambda <- log(curr_lambda) + gamma(iteration) * (alpha - alpha_opt)
    return(exp(log_lambda))
}

upd_mu <- function(curr_mu, iteration, x_next)
{
    # x_next is basically X_{i + 1}, where i is current teration.
    updated_mu <- curr_mu + gamma(iteration) * (x_next - curr_mu)
    return(updated_mu)
}

upd_diag_sigma <- function(curr_sigma, iteration, curr_mu, x_next)
{
    # note the update step is actually sigma_new = sigma_old + gamma * (outer_product(X - mu, X_mu) - sigma_old)
    # but here, we consider only diagonal sigma,so the outer product reduces to a simple element wise square.
    # we take the diagonal entries of the outer product only. This is also inline with the original
    # barkers proposal paper.(refer to page 17, first line)
    updated_sigma <- curr_sigma + gamma(iteration) * ((x_next - curr_mu)^2 - curr_sigma)
    return(updated_sigma)
}

#---------------Barker's Algo-----------------------------
barker_adap_aryt <- function(y, X, N = 1e3, dist = "normal")
{
    p <- dim(X)[2]
    alpha_opt <- 0.574

    # initalization of beta. note that glm was throwing a NA value due to small data available.
    foo <- glm(y ~ X - 1, family = binomial("logit"))$coef
    beta <- as.matrix(foo, ncol = 1)
    col_na <- which(is.na(beta))
    beta[col_na] <- colMeans(X)[col_na]
    beta.mat <- matrix(0, nrow = N, ncol = p)
    beta.mat[1, ] <- as.numeric(beta)

    #mu_inital
    mu <- rep(0, p)
    # sigma initial
    sig <- rep(1, p)
    #scale initial - refer to original barker's paper, page 17, line before section 6.2
    lambda <- (2.4 ^ 2) / (p ^ (1 / 3))

    accept <- 0

    for (i in 2:N)
    {
        grad_beta <- grad_log_posterior(beta, X, y)

        # sample diagonally!
        z <- samp_z(n = p, h = lambda * sig, dist = dist)
        prob_invert <- 1 / (1 + exp(-z * grad_beta))
        inv_or_not <- (runif(p) < prob_invert)
        b <- 2 * inv_or_not - 1
        prop <- beta + z * b

        grad_prop <- grad_log_posterior(prop, X, y)
        log_alpha <- log_posterior(prop, X, y)
        log_alpha <- log_alpha - log_posterior(beta, X, y)
        log_alpha <- - sum(log1p(exp( grad_prop * (prop - beta) ))) + sum(log1p(exp( grad_beta * (beta - prop) )))

        if (log(runif(1)) < log_alpha)
        {
            beta <- prop
            accept <- accept + 1
        }

        beta.mat[i, ] <- beta

        alpha <- min(1, exp(log_alpha))
        lambda <- upd_lambda(lambda, i, alpha, alpha_opt)
        # note, sigma update is being done before mu update, implying older value of mu is used, not newer.
        sig <- upd_diag_sigma(sig, i, mu, beta)
        mu <- upd_mu(mu, i, beta)
    }
    ret <- list(chain = beta.mat, accept = accept / N)
    return(ret)
}

sgbd_aryt <- function(y, X, N = 1e3, minibatch_size = 452, dist = "normal")
{
    p <- dim(X)[2]
    dataset_size <- nrow(X)

    # initalization of beta. note that glm was throwing a NA value due to small data available.
    foo <- glm(y ~ X - 1, family = binomial("logit"))$coef
    beta <- as.matrix(foo, ncol = 1)
    col_na <- which(is.na(beta))
    beta[col_na] <- colMeans(X)[col_na]
    beta.mat <- matrix(0, nrow = N, ncol = p)
    beta.mat[1, ] <- as.numeric(beta)

    #mu_inital
    mu <- rep(0, p)
    # sigma initial
    sig <- rep(1, p)
    #scale initial - refer to original barker's paper, page 17, line before section 6.2
    lambda <- (2.4 ^ 2) / (p ^ (1 / 3))

    for (i in 2:N)
    {
        rand_perm <- sample(dataset_size)[1:minibatch_size]
        X_mini <- X[rand_perm, ]
        y_mini <- y[rand_perm]
        z <- samp_z(n = p, h = lambda * sig, dist = dist)
        prob_invert <- 1 / (1 + exp(-z * grad_log_posterior(beta, X_mini, y_mini)))
        inv_or_not <- (runif(p) < prob_invert)
        b <- 2 * inv_or_not - 1
        beta <- beta + z * b
        beta <- beta.mat[i, ]
    }

    return(beta.mat)
}


#--------------------DATASET PREP---------------------

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