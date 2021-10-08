#' @title Update FFhat at time t
#' @description Update the variational parameter FFhat at time t
#'     in the variational distribution q(X) which is a state space model
#' @param V_inv_FF_lst list of \code{V^{-1}FF} with elements being value under different
#'     value of discrete hidden states
#' @param Vhat_inv variational parameter inverse of Vhat
#' @param qs_t probability vector of the discrete hidden state at time t
#'
#' @return An updated value of FFhat at time t
#' @export
#'

update_FFhat <- function(V_inv_FF_lst, Vhat_inv, qs_t) {
    wgt_V_inv_FF <- wt_mean(V_inv_FF_lst, qs_t)
    FFhat <- solve(Vhat_inv, wgt_V_inv_FF)
    return(FFhat)
}
