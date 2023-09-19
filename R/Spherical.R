#' @title Spherical Interpolation Function
#'
#' @description
#' Spherical interpolation (slerp) is a method for smoothly interpolating between two values on a unit sphere. 
#' This function calculates slerp for a given input vector `x` with user-defined parameters.
#'
#' @param x A numeric vector or scalar containing input values.
#' @param n The starting value on the unit sphere.
#' @param r The rotation angle in radians.
#' @param s The ending value on the unit sphere.
#'
#' @return 
#' A numeric vector of the same length as `x`, representing the slerp interpolation between `n` and `s`.
#' 
#' @export
#' 
#' @details
#' This function is adapted from the \href{https://github.com/ET-NCMP/NCMP}{ET-NCMP/NCMP} and modified for use in R-Instat.
#' The modifications include changes to data input and output processes while preserving the core calculation methods. 
#' Slerp is often used in computer graphics and animation to smoothly interpolate between two rotations or orientations. 
#' This function calculates slerp for a given input vector `x` based on the parameters `n` (start value), `r` (rotation angle), 
#' and `s` (end value). It ensures that the interpolation stays on the unit sphere.
#'
#' @examples
#' # Perform spherical interpolation between (0, 0, 1) and (1, 0, 0) with rotation of pi/2
#' # x <- seq(0, 1, by = 0.1)
#' # n <- c(0, 0, 1)
#' # r <- pi / 2
#' # s <- c(1, 0, 0)
#' # result <- Spherical(x, n, r, s)
#' # plot(result, type = "l", main = "Spherical Interpolation")
#'
#' @references 
#' For the original source code and more information, please refer to: \href{https://github.com/ET-NCMP/NCMP}{ET-NCMP/NCMP}
#' 
Spherical <- function (x, n, r, s) { ifelse(x <= r, (s - n) * (3 * x / (2 * r) - x ^ 3 / (2 * r ^ 3)) + n, s) }
