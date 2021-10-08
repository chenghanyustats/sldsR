#' @title Update What_inv at time t
#' @description Update the variational parameter What_inv at time t
#'     in the variational distribution q(X) which is a state space model
#' @param W_inv_lst list of \code{W^{-1}} with elements being value under different
#'     value of discrete hidden states
#' @param FF_V_inv_FF_lst list of FF %*% V^{-1} FF with elements being value under different
#'     value of discrete hidden states
#' @param Vhat_inv variational parameter inverse of Vhat
#' @param FFhat variational parameter FFhat
#' @param GG_W_inv_GG_lst list of GG %*% W^{-1} GG with elements being value under different
#'     value of discrete hidden states
#' @param What_inv variational parameter inverse of What
#' @param GGhat variational parameter GGhat
#' @param qs_t probability vector of the discrete hidden state at time t
#' @param t time t
#' @param TT time length
#'
#'
#' @return An updated value of What inverse at time t
#' @export
#'

update_What_inv <- function(W_inv_lst,
                            FF_V_inv_FF_lst, Vhat_inv, FFhat,
                            GG_W_inv_GG_lst = NULL,
                            What_inv = NULL,
                            GGhat = NULL,
                            qs_t, t, TT) {
    A <- wt_mean(W_inv_lst, qs_t) +
        wt_mean(FF_V_inv_FF_lst, qs_t) -
        emulator::quad.form(Vhat_inv, FFhat)
    if (t == TT) {
        return(A)
    } else {
        B <- A + wt_mean(GG_W_inv_GG_lst, qs_t) -
            emulator::quad.form(What_inv, GGhat)
        return(B)
    }
}
