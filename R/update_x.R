#' @title Update continous latent states using FFBS algorithm
#' @description Sample new continous latent states using forward filtering
#'     backward sampling (FFBS) algorithm with the normal prior \code{N(m0, C0)}
#' @param m0 mean of the normal prior
#' @param C0 (co)variance of the normal prior
#' @param y observation vector
#' @param v mcmc sample of observational variance
#' @param w mcmc sample of AR(1) coefficient
#' @param phi mcmc sample of AR(1) coefficient
#' @param st mcmc sample of \code{s_t}, the discrete hidden states
#' @param k2 observational variance when \code{s} is 2
#'
#' @return A vector of continous latent state
#' @export
#'



update_x <- function(m0, C0, y, v, w, phi, st, k2) {
    k2_idx <- st
    k2_idx[st == 2] <- k2
    mod <- dlm::dlm(FF = 1, V = v, JV = 1,
               GG = phi, W = w,
               m0 = m0, C0 = C0,
               X = k2_idx * v)
    filt <- dlm::dlmFilter(y, mod)
    x <- dlm::dlmBSample(filt)
    # print(filt)
    # print(x)
    # while(Inf %in% x || -Inf %in% x) {
    #     x <- dlm::dlmBSample(filt)
    # }
    return(as.vector(x))
}
