#' @title log forward probabilities of HMM
#' @description Compute (log) forward probability of a hidden Markov model
#' @param K the number of switching states
#' @param n the number of time points
#' @param h_mat likelihood value matrix
#' @param init_dist initial dostribution
#' @param tran_prob transition probability matrix
#'
#' @return A log forward probability vector
#' @export
#'

log_forward_hmm <- function(K, n, h_mat, init_dist, tran_prob){
    log_alpha <- matrix(0, K, n)
    foo <- init_dist * h_mat[, 1]
    w <- sum(foo)
    ll <- log(w)
    phi <- foo / w
    log_alpha[, 1] <- ll + log(foo)

    for(i in 2:n) {
        foo <- phi %*% tran_prob * h_mat[, i]
        w <- sum(foo)
        ll <- ll + log(w)
        phi <- foo / w
        log_alpha[, i] <- ll + log(phi)
    }
    return(log_alpha)
}
