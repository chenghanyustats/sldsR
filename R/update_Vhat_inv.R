#' @title Update inverse of Vhat at time t
#' @description Update the variational parameter Vhat inverse at time t
#'     in the variational distribution q(X) which is a state space model
#' @param V_inv_lst list of inverse of V with elements being value under different
#'     value of discrete hidden states
#' @param qs_t probability vector of the discrete hidden state at time t
#'
#' @return An updated value of inverse of Vhat at time t
#' @export
#'

update_Vhat_inv <- function(V_inv_lst, qs_t) {
    K <- length(qs_t)
    summ <- 0
    for (i in 1:K) {
        summ <- summ + V_inv_lst[[i]] * qs_t[i]
    }
    return(summ)
}
