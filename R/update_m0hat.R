#' @title Update m0hat
#' @description Update the variational parameter m0hat
#'     in the variational distribution q(X) which is a state space model
#' @param C0hat_inv inverse of C0hat
#' @param C0_inv_m0_lst list of \code{C_0^{-1}m_0} with elements being value under different
#'     value of discrete hidden states
#' @param qs_t probability vector of the discrete hidden state at time t
#'
#'
#' @return An updated value of m0hat
#' @export
#'

update_m0hat <- function(C0hat_inv, C0_inv_m0_lst, qs_t) {
    wgt_C0_inv_m0 <- wt_mean(C0_inv_m0_lst, qs_t)
    return(solve(C0hat_inv, wgt_C0_inv_m0))
}
