# Math required for pulse-specific, weighted relative amplitude
# estimation.

library(stats)


#' Returns a penalization factor for some z-score, such that
#'
#' @param x Subject ID
quad_penalty <- function(z, alpha = 1) {
    return(alpha / (z^2 + 1))
}

gaussian_penalty <- function(z, alpha = 1) {
    return(exp(-z^2 / alpha))
}
gaussian_penalty(-2, alpha = 10)
quad_penalty(3, alpha = 4)

robust_zscores <- function(v) {
    mad <- mad(v)
    rz <- 0.6745 * (v - median(v)) / mad

    return(rz)
}

