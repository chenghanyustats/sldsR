#' @title Update GGhat at time t
#' @description Update the variational parameter GGhat at time t
#'     in the variational distribution q(X) which is a state space model
#' @param What_inv variational parameter inverse of What
#' @param W_inv_GG_lst list of \code{W^{-1} GG} with elements being value under different
#'     value of discrete hidden states
#' @param qs_t probability vector of the discrete hidden state at time t
#'
#' @return An updated value of GGhat at time t
#' @export
#'



update_GGhat <- function(What_inv, W_inv_GG_lst, qs_t) {
    wgt_W_inv_GG <- wt_mean(W_inv_GG_lst, qs_t)
    GGhat <- solve(What_inv, wgt_W_inv_GG)
}
