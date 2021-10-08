#' @title Update AR(1) coefficient in the evolution level process
#' @description Sample AR(1) coefficient from a normal distribution
#'     with a flat prior with a constriant that the coefficient is in \code{[-1, 1]}.
#' @param xt mcmc sample of \code{x_t}, the continuous latent states
#' @param xm1 mcmc sample of \code{x_{t-1}}, the continuous latent states
#' @param w mcmc sample of evolutional variance
#'
#' @return A number drawn from a normal distribtion
#' @export
#'


update_phi <- function(xt, xm1, w) {
    sum_Xm1_2 <- crossprod(xm1)
    m_phi <- crossprod(xt, xm1) / sum_Xm1_2
    C_phi <- w / sum_Xm1_2
    phi <- stats::rnorm(1, mean = m_phi, sd = sqrt(C_phi))
    while (phi > 1 || phi < -1) {
        phi <- stats::rnorm(1, mean = m_phi, sd = sqrt(C_phi))
    }
    return(phi)
}
