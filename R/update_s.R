#' @title Update discrete hidden states using discrete version of FFBS algorithm
#' @description Sample new discrete hidden states using discrete version of
#'     forward filtering backward sampling (FFBS) algorithm with a given transition
#'     matrix.
#' @param y observation vector
#' @param xt mcmc sample of \code{x_t}, the continuous latent states
#' @param v mcmc sample of observational variance
#' @param k2 observational variance when \code{s} is 2
#' @param tran_mat transition matrix
#' @return A vector of discrete hidden state
#' @export
#'

update_s <- function(y, xt, v, k2, tran_mat) {
    TT <- length(y)
    p_ahead_1 <- rep(0, TT)
    p_filter_1 <- rep(0, TT)
    p_0_1 <- 0.5
    p_0_k2 <- 1 - p_0_1

    sd_v <- sqrt(v)
    sd_v_k2 <- sqrt(k2 * v)

    ## First one
    # 1. One-step ahead
    # =================
    p_ahead_1[1] <- tran_mat[1, 1] * p_0_1 + tran_mat[2, 1] * p_0_k2

    # 2. filtering
    # =================
    A1 <- stats::dnorm(y[1], xt[1], sd_v, log = TRUE) + log(p_ahead_1[1])
    A2 <- stats::dnorm(y[1], xt[1], sd_v_k2, log = TRUE) + log(1 - p_ahead_1[1])
    p_filter_1[1] <- exp(A1) / (exp(A1) + exp(A2))

    ## Second to the last
    for (t in 2:TT) {
        # 1. One-step ahead
        # =================
        p_ahead_1[t] <- tran_mat[1, 1] * p_filter_1[t - 1] +
            tran_mat[2, 1] * (1 - p_filter_1[t - 1])

        # 2. filtering
        # =================
        A1 <- stats::dnorm(y[t], xt[t], sd_v, log = TRUE) + log(p_ahead_1[t])
        A2 <- stats::dnorm(y[t], xt[t], sd_v_k2, log = TRUE) + log(1 - p_ahead_1[t])
        p_filter_1[t] <- exp(A1) / (exp(A1) + exp(A2))
    }

    St_vec <- rep(0, TT)
    ###################
    # sample S_T
    ###################
    ## Last one
    rT <- stats::rmultinom(1, size = 1, prob = c(p_filter_1[TT], 1 - p_filter_1[TT]))
    St_vec[TT] <- rT[1] + 2 * rT[2]
    ## Second last to the first
    for (t in (TT - 1):1) {
        if(St_vec[t + 1] == 1) {
            B1 <- tran_mat[1, 1] * p_filter_1[t]
            B2 <- tran_mat[2, 1] * (1 - p_filter_1[t])
            rt <- stats::rmultinom(1, size = 1, prob = c(B1, B2) / p_ahead_1[t + 1])
            St_vec[t] <- rt[1] + 2 * rt[2]
        } else {
            B1 <- tran_mat[1, 2] * p_filter_1[t]
            B2 <- tran_mat[2, 2] * (1 - p_filter_1[t])
            rt <- stats::rmultinom(1, size = 1, prob = c(B1, B2) / (1 - p_ahead_1[t + 1]))
            St_vec[t] <- rt[1] + 2 * rt[2]
        }
    }
    return(St_vec)
}
