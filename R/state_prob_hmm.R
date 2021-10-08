#' @title posterior state probabilities of HMM
#' @description Compute posterior discrete state probabilities of a hidden Markov model
#' @param K the number of switching states
#' @param n the number of time points
#' @param h_mat likelihood value matrix
#' @param init_dist initial dostribution
#' @param tran_prob transition probability matrix
#'
#' @return A posterior state probability vector
#' @export
#'

state_prob_hmm <- function(K, n, h_mat, init_dist, tran_prob) {
    log_alpha <- log_forward_hmm(K, n, h_mat, init_dist, tran_prob)
    log_beta <- log_backward_hmm(K, n, h_mat, init_dist, tran_prob)
    c <- max(log_alpha[, n])
    log_lik <- c + log(sum(exp(log_alpha[, n] - c)))
    state_prob <- matrix(0, nrow = K, ncol = n)
    for (i in 1:n) {
        state_prob[, i] <- exp(log_alpha[, i] + log_beta[, i] -
                                   log_lik)
    }
    return(state_prob)
}
