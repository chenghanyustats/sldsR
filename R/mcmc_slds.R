#' @title Main MCMC algorithm for SLDS
#' @description Generate posterior samples of unknown continuous and discrete latent
#'     states as well as other unknown model parameters.
#' @param y observation vector
#' @param x_init initial value of continuous latent states
#' @param s_init initial value of discrete hidden states
#' @param v_init initial value of observational variance
#' @param w_init initial value of evolution variance
#' @param phi_init initial value of AR(1) coefficient
#' @param xi_init initial probability of the transition matrix
#' @param av0 shape parameter of the inverse gamma prior for v
#' @param bv0 scale parameter of the inverse gamma prior for v
#' @param aw0 shape parameter of the inverse gamma prior for w
#' @param bw0 scale parameter of the inverse gamma prior for w
#' @param a_xi first shape parameter of the beta prior for xi
#' @param b_xi second shape parameter of the beta prior for xi
#' @param m0 mean of the normal prior
#' @param C0 (co)variance of the normal prior
#' @param k2 observational variance when \code{s} is 2
#' @param burn burning period
#' @param thin thinning
#' @param mcmc_len total number of MCMC iterations
#'
#' @return A list of posterior samples of of unknown continuous and discrete latent
#'     states and unknown model parameters
#' @export
#'

mcmc_slds <- function(y,
                      x_init = stats::rnorm(length(y) + 1),
                      s_init = sample(c(1, 2), size = length(y) + 1, replace = TRUE),
                      v_init = 2,
                      w_init = 2,
                      phi_init = 0.5,
                      xi_init = 0.5,
                      av0 = 1 / 2, bv0 = 1 / 2,
                      aw0 = 1 / 2, bw0 = 1 / 2,
                      a_xi = 1, b_xi = 1,
                      m0 = 0, C0 = 1000,
                      k2 = 30,
                      burn = 0,
                      thin = 1,
                      mcmc_len = 2000) {

    TT <- length(y)

    # =========================================================
    ## storage of samples of parameters
    # =========================================================
    # idx recording which iterations are sampled
    sampleidx <- seq(from = (burn + thin), to = mcmc_len, by = thin)
    # sample storage
    sample_len <- length(sampleidx)

    V <- rep(NA, sample_len)
    W <- rep(NA, sample_len)
    Phi <- rep(NA, sample_len)
    Xi <- rep(NA, sample_len)
    X <- matrix(NA, nrow = TT + 1, ncol = sample_len)
    St <-  matrix(NA, nrow = TT, ncol = sample_len)


    # =========================================================
    ## initial values
    # =========================================================
    # X[, 1] <- rnorm(TT + 1)
    # Xt <- X[-1, ]
    # Xm1 <- X[-(TT + 1), ]
    # a <- rmultinom(1, size = 1, prob = c(0.5, 0.5))
    # S[1, 1] <- a[1] + 2 * a[2]
    # for (i in 2:(TT + 1)) {
    #     if(S[i - 1, 1] == 1) {
    #         r1 <- rmultinom(1, size = 1, prob = c(p, 1 - p))
    #         S[i, 1] <- r1[1] + 2 * r1[2]
    #     } else {
    #         r2 <- rmultinom(1, size = 1, prob = c(1 - p, p))
    #         S[i, 1] <- r2[1] + 2 * r2[2]
    #     }
    # }
    # St <- S[-1, ]
    # V[1] <- 2
    # W[1] <- 2
    # Phi[1] <- 0.5

    x <- x_init
    xt <- x[-1]
    xm1 <- x[-(TT + 1)]
    s <- s_init
    st <- s[-1]
    v <- v_init
    w <- w_init
    phi <- phi_init
    xi <- xi_init
    tran_mat <- matrix(c(xi, 1 - xi, 1 - xi, xi), 2, 2)

    # =========================================================
    ## MCMC iteration
    # =========================================================
    for (mcmc_iter in 1:mcmc_len) {
        if (mcmc_iter %% 100 == 0) cat("MCMC iter:", mcmc_iter, "\r")

        #=========
        # sample v
        #=========
        # av <- av0 + TT / 2
        # idx_1 <- which(st == 1)
        # idx_k2 <- which(st == 2)
        # bv <- bv0 + crossprod(y[idx_1] - xt[idx_1]) / 2 +
        #     crossprod(y[idx_k2] - xt[idx_k2]) / (2 * k2)
        # (v <- 1 / rgamma(1, shape = av, rate = bv))


        v <- update_v(av0 = av0, bv0 = bv0, y = y, xt = xt, st = st, k2 = k2)
        # print(paste("v=", v))
        #=========
        # sample w
        #=========
        w <- update_w(aw0 = aw0, bw0 = bw0, xt, xm1, phi)
        # print(paste("w=", w))

        #=========
        # sample phi
        #=========
        phi <- update_phi(xt = xt, xm1 = xm1, w = w)
        # print(paste("phi=", phi))

        #=========
        # sample x using FFBS (dlm package)
        #=========
        # for(i in 1:100) {
        #     print(x <- update_x(m0 = m0, C0 = C0, y = y, v = v, w = w, phi = phi,
        #                    st = st, k2 = k2))
        # }
        x <- update_x(m0 = m0, C0 = C0, y = y, v = v, w = w, phi = phi,
                      st = st, k2 = k2)
        xt <- as.vector(x[-1])
        xm1 <- as.vector(x[-(TT + 1)])
        # print(paste("x=", x))
        #=========
        # sample S
        #=========
        st <- update_s(y = y, xt = xt, v = v, k2 = k2, tran_mat = tran_mat)
        # st <- s[-1]
        # print(paste("s=", s))
        #=========
        # sample xi and update transition matrix
        #=========
        xi_res <- update_xi(a_xi = a_xi, b_xi = b_xi, st = st)
        xi <- xi_res$xi
        tran_mat <- xi_res$tran_mat

        #=========
        # save samples
        #=========
        # V[mcmc_iter] <- v
        # W[mcmc_iter] <- w
        # Phi[mcmc_iter] <- phi
        # Xi[mcmc_iter] <- xi
        # X[, mcmc_iter] <- x
        # S[, mcmc_iter] <- s
        iter <- mcmc_iter - burn
        if (mcmc_iter > burn) {
            iter <- mcmc_iter - burn
            if (mcmc_iter %% thin == 0) {
                V[iter %/% thin] <- v
                W[iter %/% thin] <- w
                Phi[iter %/% thin] <- phi
                Xi[iter %/% thin] <- xi
                X[, iter %/% thin] <- x
                St[, iter %/% thin] <- st
            }
        }
    }
    return(list(V = V,
                W = W,
                Phi = Phi,
                Xi = Xi,
                X = X,
                St = St,
                burn = burn,
                thin = thin,
                mcmc_len = mcmc_len,
                sampleidx = sampleidx))
}
