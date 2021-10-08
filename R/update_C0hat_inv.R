#' @title Update C0hat_inv
#' @description Update the variational parameter C0hat_inv
#'     in the variational distribution q(X) which is a state space model
#' @param C0_inv_lst list of \code{C_0^{-1}} with elements being value under different
#'     value of discrete hidden states
#' @param FF_V_inv_FF_lst list of FF' V^-1 FF with elements being value under different
#'     value of discrete hidden states
#' @param Vhat_inv variational parameter inverse of Vhat
#' @param FFhat variational parameter FFhat
#' @param GG_W_inv_GG_lst list of GG' W^-1 GG with elements being value under different
#'     value of discrete hidden states
#' @param What_inv variational parameter inverse of What
#' @param GGhat variational parameter GGhat
#' @param qs_t probability vector of the discrete hidden state at time t
#'
#'
#' @return An updated value of C0hat inverse
#' @export
#'


update_C0hat_inv <- function(C0_inv_lst, FF_V_inv_FF_lst,
                             Vhat_inv, FFhat,
                             GG_W_inv_GG_lst,
                             What_inv,
                             GGhat,
                             qs_t) {
    C0hat_inv <- wt_mean(C0_inv_lst, qs_t) +
        wt_mean(FF_V_inv_FF_lst, qs_t) -
        emulator::quad.form(Vhat_inv, FFhat) +
        wt_mean(GG_W_inv_GG_lst, qs_t) -
        emulator::quad.form(What_inv, GGhat)
    return(C0hat_inv)
}
