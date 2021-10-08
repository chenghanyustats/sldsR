################################################################################
# Implementation of SLDS and Analysis                                          #
# Cheng-Han Yu                                                                 #
################################################################################
rm(list = ls())
################################################################################
## Load package and data
################################################################################
library(sldsR)
load("./data/ar1_noise_mixture_data.rda", verbose = TRUE)
################################################################################
## Run the MCMC algorithm
################################################################################
y <- ar1_noise_mixture_data$y
TT <- ar1_noise_mixture_data$TT
k2 <- ar1_noise_mixture_data$k2

system.time(mcmc_result <- mcmc_slds(y = y,
                         x_init = stats::rnorm(TT + 1),
                         s_init = sample(c(1, 2), size = TT + 1, replace = TRUE),
                         k2 = 30,
                         burn = 500, thin = 1, mcmc_len = 2000))

################################################################################
## Posterior inference
################################################################################
### trace plots and histograms
par(mfrow = c(4, 2))
par(mar = c(4, 4, 2, 1))
plot(mcmc_result$Phi, type = 'l')
abline(h = ar1_noise_mixture_data$phi, col = 2)
hist(mcmc_result$Phi, breaks = 30, xlab = "phi1", main = "Histogram of phi")
plot(mcmc_result$V, type = 'l')
abline(h = ar1_noise_mixture_data$v, col = 2)
hist(mcmc_result$V, breaks = 30, xlab = "V", main = "Histogram of v")
plot(mcmc_result$W, type = 'l')
abline(h = ar1_noise_mixture_data$w, col = 2)
hist(mcmc_result$W, breaks = 30, xlab = "W", main = "Histogram of w")
plot(mcmc_result$Xi, type = 'l')
abline(h = ar1_noise_mixture_data$xi, col = 2)
hist(mcmc_result$Xi, breaks = 30, xlab = "W", main = "Histogram of xi")


### filtering and posterior probabilty
par(mfrow = c(2, 1))
par(mar = c(4, 4, 4, 1))
plot(y, xlab = "time", axes = FALSE, pch = 19, col = "green3",
     main = expression(paste("data (circles) and filtering p(", x[t], " | ", y[1:t], ")")))
axis(1)
axis(2)
idx <- which(ar1_noise_mixture_data$s[2:TT] == k2)
points(idx, y[idx], pch = 19)

posmean <- apply(mcmc_result$X, 1, mean)
lower <- apply(mcmc_result$X, 1, quantile, prob = 0.025)
upper <- apply(mcmc_result$X, 1, quantile, prob = 0.975)
lines(posmean, col = 2, lwd = 2)
lines(lower, col = "blue", lwd = 1.5, lty = 1)
lines(upper, col = "blue", lwd = 1.5, lty = 1)

prob_st <- apply(mcmc_result$St, 1, function(x) {
    sum(x == 1) / length(mcmc_result$sampleidx)
})
plot(1:TT, prob_st, axes = FALSE, pch = 20, xlab = "time", col = "green3",
     ylab = expression(paste("Pr(", S[t], " = 1 | ", y[1:200], ")")),
     main = expression(paste("Pr(", S[t], " = 1 | ", y[1:200], ")")))
points(idx, prob_st[idx], pch = 20)
abline(h = 0.5, col = "red", lty = 2)
axis(1)
axis(2)

## posterior distributions and 95% intervals
par(mfrow = c(2, 2))
par(mar = c(4, 4, 4, 1))
hist(mcmc_result$Phi, prob = F, breaks = 30, border = FALSE, col = "gray",
     main = expression(paste("Posterior of ", phi)), xlab = "")
abline(v = ar1_noise_mixture_data$phi, col = 2, lwd = 2)
abline(v = quantile(mcmc_result$Phi, prob = 0.025), col = 4, lwd = 2)
abline(v = quantile(mcmc_result$Phi, prob = 0.975), col = 4, lwd = 2)

hist(mcmc_result$V, prob = F, breaks = 30, border = FALSE, col = "gray",
     main = expression(paste("Posterior of v")), xlab = "")
abline(v = ar1_noise_mixture_data$v, col = 2, lwd = 2)
abline(v = quantile(mcmc_result$V, prob = 0.025), col = 4, lwd = 2)
abline(v = quantile(mcmc_result$V, prob = 0.975), col = 4, lwd = 2)

hist(mcmc_result$W, prob = F, breaks = 30, border = FALSE, col = "gray",
     main = expression(paste("Posterior of w")), xlab = "")
abline(v = ar1_noise_mixture_data$w, col = 2, lwd = 2)
abline(v = quantile(mcmc_result$W, prob = 0.025), col = 4, lwd = 2)
abline(v = quantile(mcmc_result$W, prob = 0.975), col = 4, lwd = 2)

hist(mcmc_result$Xi, prob = F, breaks = 30, border = FALSE, col = "gray",
     main = expression(paste("Posterior of ", xi)), xlab = "")
abline(v = ar1_noise_mixture_data$xi, col = 2, lwd = 2)
abline(v = quantile(mcmc_result$Xi, prob = 0.025), col = 4, lwd = 2)
abline(v = quantile(mcmc_result$Xi, prob = 0.975), col = 4, lwd = 2)



################################################################################
## Run the VI algorithm
################################################################################
vi_result <- variational_em_slds(y = y, qs = matrix(0.5, 2, length(y)))

################################################################################
## Posterior inference
################################################################################
par(mfrow = c(2, 1))
par(mar = c(4, 4, 4, 1))
plot(y, xlab = "time", axes = FALSE, pch = 19, col = "green3",cex = 0.8,
     main = expression(paste("data (circles) and filtering p(", x[t], " | ", y[1:t], ") Iter = 100")))
axis(1)
axis(2)
# idx <- which(St_vec[2:TT] == k2)
points(idx, y[idx], pch = 19,cex = 0.8)

lines(1:TT, vi_result$ms, col = 2, lwd = 2)
lines(1:TT, vi_result$ms + 1.96 * sqrt(as.vector(vi_result$Cs)), col = "blue", lwd = 1.5, lty = 1)
lines(1:TT, vi_result$ms - 1.96 * sqrt(as.vector(vi_result$Cs)), col = "blue", lwd = 1.5, lty = 1)

plot(1:TT, vi_result$qs[2, ], axes = FALSE, pch = 20, xlab = "time", col = "green3",
     ylab = expression(paste("Pr(", S[t], " = 1 | ", y[1:TT], ")")),
     main = expression(paste("Pr(", S[t], " = 1 | ", y[1:TT], ")")))
points(idx, vi_result$qs[2, idx], pch = 20)
abline(h = 0.5, col = "red", lty = 2)
axis(1)
axis(2)

par(mfrow = c(1, 1))
plot(vi_result$ELBO_vec, type = "l", main = "Evidence Lower Bound", xlab = "Iterations")
