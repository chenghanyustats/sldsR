#' @title Compute weighted mean
#' @description Compute weighted means from a given list and weight vector
#' @param x_lst a list containing elements to be weighted
#' @param wt weight vector
#'
#' @return A weigthed mean of a given list
#' @export
#'


wt_mean <- function(x_lst, wt) {
    K <- length(wt)
    summ <- 0
    for (i in 1:K) {
        summ <- summ + x_lst[[i]] * wt[i]
    }
    return(summ)
}
