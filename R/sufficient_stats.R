#' @title sufficient statistics compuation
#' @description Compute sufficient statistics of a state space model
#' @param ms mean of Kalman smoother
#' @param Cs (co)variance of Kalman smoother
#' @param tm1 time t - 1
#' @param t time t
#' @param GG evolution matrix
#' @param m0 mean of initial continuous latent state
#' @param C0 mean of initial continuous latent state
#' @return A sufficient statistics value
#' @export
#'




sufficient_stats <- function(ms, Cs, tm1, t, GG = NULL, m0 = NULL, C0 = NULL) {
    if (tm1 == t) {
        # t <- t1
        Exx <- Cs[, , t] + ms[, t] %*% t(ms[, t])
        return(Exx)
    } else {
        # if (tm1 == 0) {
        #     Ex1x <- (C0 + m0 %*% t(m0)) %*% t(GG) + m0 %*% t(ms[, t])
        # } else {
        #     Ex1x <- (Cs[, , tm1] + ms[, tm1] %*% t(ms[, tm1])) %*% t(GG) + ms[, tm1] %*% t(ms[, t])
        # }
        if (tm1 == 0) {
            Ex1x <- C0 %*% t(GG) + m0 %*% t(ms[, t])
        } else {
            Ex1x <- Cs[, , tm1] %*% t(GG) + ms[, tm1] %*% t(ms[, t])
        }
        return(Ex1x)
    }
}
