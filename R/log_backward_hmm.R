#' @title log backward probabilities of HMM
#' @description Compute (log) backward probability of a hidden Markov model
#' @param K the number of switching states
#' @param n the number of time points
#' @param h_mat likelihood value matrix
#' @param init_dist initial dostribution
#' @param tran_prob transition probability matrix
#'
#' @return A log backward probability vector
#' @export
#'

log_backward_hmm <- function(K, n, h_mat, init_dist, tran_prob){
    log_beta <- matrix(0, K, n)
    log_beta[, n] <- rep(0, K)
    foo <- rep(1/K, K)
    ll <- log(K)
    for (i in (n - 1):1) {
        foo <- tran_prob %*% (h_mat[, i + 1] * foo)
        log_beta[, i] <- log(foo) + ll
        w <- sum(foo)
        phi <- foo / w
        ll <- ll + log(phi)
    }
    return(log_beta)
}
