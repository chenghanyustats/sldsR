#' @title Update probability in the transition matrix
#' @description Sample probability from a beta distribution
#'     with a \code{beta(a_xi, b_xi)} prior
#' @param a_xi first shape parameter of the beta prior
#' @param b_xi second shape parameter of the beta prior
#' @param st mcmc sample of \code{s_t}, the discrete hidden states
#'
#' @return A 2 by 2 transition matrix
#' @export
#'

update_xi <- function(a_xi, b_xi, st) {
    # N_stay <- 0
    # N_vary <- 0
    # for (i in 2:length(st)) {
    #     if (st[i] == st[i - 1]) {
    #         N_stay <- N_stay + 1
    #     } else {
    #         N_vary <- N_vary + 1
    #     }
    # }
    #
    N_stay <- sum(st[-1] == st[-length(st)])
    N_vary <- length(st) - N_stay - 1
    xi <- stats::rbeta(1, shape1 = N_stay + a_xi, shape2 = N_vary + b_xi)
    return(list(xi = xi,
                tran_mat = matrix(c(xi, 1 - xi, 1 - xi, xi), 2, 2)))
}
