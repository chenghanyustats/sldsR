## code to prepare `ar1_noise_mixture` dataset goes here
set.seed(12345)
TT <- 500
phi <- 0.9
k2 <- 30
v <- 1
w_var <- 0.5
p <- 0.99
tran_mat <- matrix(c(p, 1 - p, 1 - p, p), 2, 2)
a_true <- rmultinom(1, size = 1, prob = c(0.5, 0.5))
S0 <- v * (a_true[1] + k2 * a_true[2])
St_vec <- numeric(TT + 1)
St_vec[1] <- S0
for (i in 2:length(St_vec)) {
    if(St_vec[i - 1] == 1) {
        r1 <- rmultinom(1, size = 1, prob = c(p, 1 - p))
        St_vec[i] <- v * (r1[1] + k2 * r1[2])
    } else {
        r2 <- rmultinom(1, size = 1, prob = c(1 - p, p))
        St_vec[i] <- v * (r2[1] + k2 * r2[2])
    }
}

e <- rnorm(TT, 0, sqrt(St_vec[2:TT]))
w <- rnorm(TT, 0, sqrt(w_var))
y <- numeric(TT)
x <- numeric(TT)
x0 <- 0 # mu_0 = 0
x[1] <- phi * x0 + w[1]
y[1] <- x[1] + e[1]
for(i in 2:TT){
    x[i] <- phi * x[i - 1] + w[i]
    y[i] <- x[i] + e[i]
}

ar1_noise_mixture_data <- list(x = x,
                               y = y,
                               s = St_vec,
                               TT = TT,
                               phi = phi,
                               k2 = k2,
                               xi = p,
                               w = w_var,
                               v = v,
                               tran_mat = tran_mat)

### plotting data
# par(mfrow = c(1, 1))
# par(mar = c(4, 4, 1, 1))
# plot(y, xlab = "time", axes = FALSE)
# axis(1)
# axis(2)
# idx <- which(St_vec == k2)
# points(idx, y[idx], pch = 19)

usethis::use_data(ar1_noise_mixture_data, compress = "xz")
