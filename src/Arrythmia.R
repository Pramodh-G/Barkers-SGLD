source("utils/arrhythmia_utils.R")

# Note that we already define X, and y in the corresponding utils section.
set.seed(42)
nplot <- sample(1:p, 6)

N <- 1e3
# for functions involving prior_sig, please take a look at arrhytmia_utils.R
# This was set to 5 at the behest of the authors of the paper.
prior_sig <- 5

# Barker methods
pilot_result <- barker_aryt(y, X, N = 1e4, h = 1.05e-1, dist = "normal")
accept_barker <- pilot_result$accept
print(accept_barker)
pre_sd <- sqrt(diag(cov(pilot_result$chain)))
print(pre_sd)

# result <- barker_aryt(y, X, N = N, h = 2e-1 * pre_sd, dist = "normal")
# chain_barker <- result$chain
# accept_barker <- result$accept
# print(accept_barker)

common_h <- 1e-1
chain_sgbd <- sgbd_aryt(y, X, N = N, h = common_h * pre_sd, dist = "normal")
chain_sgbd_grad <- prep_gradients(chain_sgbd, X, y)
ksd(chain_sgbd, chain_sgbd_grad)

# result <- barker_adap_aryt(y, X, N = N, dist = "normal")
# plot(result$lambdas)
# chain_ad_barker <- result$chain
# accept_ad_barker <- result$accept
# mu_ad_barker <- result$mu
# sig_ad_barker <- result$sig
# lambda_ad_barker <- result$lambda
# print(accept_ad_barker)
# print(mu_ad_barker)
# print(sig_ad_barker)
# print(lambda_ad_barker)

# Langevin methods
# pilot <- sgld_aryt(y, X, N = N, h = 0.5, dist = "normal")
# pre_sd <- sqrt(diag(cov(pilot)))
chain_sgld <- sgld_aryt(y, X, N = N, h = common_h * pre_sd, dist = "normal")
chain_sgld_grad <- prep_gradients(chain_sgld, X, y)
ksd(chain_sgld, chain_sgld_grad)
# sum(which(is.na(chain_sgld)))

# result <- mala_aryt(y, X, N = 1e4, h = 0.1 * pre_sd)
# chain_mala <- result$chain
# accept_mala <- result$accept
# print(accept_mala)

# MH Methods
# pilot_result <- mh_aryt(y, X, N = N, h = 7.4e-2, dist = "normal")
# accept_mh <- pilot_result$accept
# print(accept_mh)
# pre_sd <- sqrt(diag(cov(pilot_result$chain)))
# result <- mh_aryt(y, X, N = N, h = 1.5e-1 * pre_sd, dist = "normal")
# chain_mh <- result$chain
# accept_mh <- result$accept
# print(accept_mh)

# Load presaved MH
chain_ad_barker_load <- readRDS("variables/aryt-barker.rds")

#choose dimensions at random
# plot.ts(chain[, 1:nplot])

# pdf("plots/aryt_h_3e0.pdf")
par(mfrow = c(3, 2))
# for (i in 1:nplot)
# {
#     acf(chain[, i], main = paste("ACF for component ", i))
# }
for (i in nplot)
{
    plot(density(chain_ad_barker_load[, i]), col = "black", main = paste("density for component ", i))
    lines(density(chain_sgld[, i]), col = "green")
    lines(density(chain_sgbd[, i]), col = "red")
    legend(x = "topleft", c("Barker's","SGLD","SGBD"),
        fill = c("black","green","red"), title = "Legend")
}

# dev.off()
