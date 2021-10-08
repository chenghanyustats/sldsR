#' @title Update variance of observational noise
#' @description Sample new observational variance from inverse gamma distribution
#'     with the inverse gamma prior \code{IG(av0, bv0)}.
#' @param av0 shape parameter of the inverse gamma prior
#' @param bv0 scale parameter of the inverse gamma prior
#' @param y observation vector
#' @param xt mcmc sample of \code{x_t}, the continuous latent states
#' @param st mcmc sample of \code{s_t}, the discrete hidden states
#' @param k2 observational variance when \code{s} is 2
#'
#' @return A number drawn from an inverse gamma distribtion
#' @export
#'


update_v <- function(av0, bv0, y, xt, st, k2) {
    # av <- av0 + length(y) / 2
    idx <- st == 1
    # idx_1 <- which(idx)
    # idx_k2 <- which(!idx)
    # bv <- bv0 + crossprod(y[idx_1] - xt[idx_1]) / 2 +
    #     crossprod(y[idx_k2] - xt[idx_k2]) / (2 * k2)
    bv <- bv0 + (crossprod(y[idx] - xt[idx]) +
        crossprod(y[!idx] - xt[!idx]) / k2) / 2
    # return(1 / stats::rgamma(1, shape = av, rate = bv))
    return(1 / stats::rgamma(1, shape = av0 + length(y) / 2, rate = bv))
}
