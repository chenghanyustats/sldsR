#' @title AR (1) with normal mixture on observational noises
#'
#' @description A list of variables generated from a swiching linear dynamical system
#' with observational process \eqn{y_t = x_t + e_t}, \eqn{e_t ~ N(0, 1)} if
#' \eqn{S_t = 1} and \eqn{e_t ~ N(0, 20)} if \eqn{S_t = 2}. The evolution process is
#' \eqn{x_t = 0.9x_{t-1} + w_t}, \eqn{w ~ N(0, 1)}. The transition matrix is a symmetric
#' 2 by 2 matrix with diagonal element 0.95.
#'
#'
#' @format A list with the following elements, which are:
#' \describe{
#' \item{x}{continous latent states}
#' \item{y}{observed outputs}
#' \item{s}{discrete hidden states}
#' \item{TT}{time length}
#' \item{phi}{AR coefficient of the evolution process}
#' \item{k2}{noise variance value when the hidden state is 2}
#' \item{tran_mat}{transition matrix}
#' }
"ar1_noise_mixture_data"
