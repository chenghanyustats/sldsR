#' @title Update variance of evolution noise
#' @description Sample new evolutional variance from inverse gamma distribution
#'     with the inverse gamma prior \code{IG(aw0, bw0)}.
#' @param aw0 shape parameter of the inverse gamma prior
#' @param bw0 scale parameter of the inverse gamma prior
#' @param xt mcmc sample of \code{x_t}, the continuous latent states
#' @param xm1 mcmc sample of \code{x_{t-1}}, the continuous latent states
#' @param phi mcmc sample of AR(1) coefficient
#'
#' @return A number drawn from an inverse gamma distribtion
#' @export
#'


update_w <- function(aw0, bw0, xt, xm1, phi) {
    # aw <- aw0 + length(xt) / 2
    # bw <- bw0 + crossprod(xt - phi * xm1) / 2
    # return(1 / stats::rgamma(1, shape = aw, rate = bw))
    return(1 / stats::rgamma(1, shape = aw0 + length(xt) / 2,
                             rate = bw0 + crossprod(xt - phi * xm1) / 2))
}
