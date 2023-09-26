#' @title Exponential Distribution Calculator
#' 
#' @description
#' Calculate the probability values of an exponential distribution with user-defined parameters. 
#' 
#' @param x A numeric vector or scalar containing input values.
#' @param n The minimum value of the distribution.
#' @param r The rate parameter, controlling the rate of decay.
#' @param s The scale factor, affecting the overall shape.
#'
#' @return 
#' A numeric vector of the same length as `x`, representing the probability values of the exponential distribution.
#' 
#' @export
#'
#' @details
#' This function is adapted from the \href{https://github.com/ET-NCMP/NCMP}{ET-NCMP/NCMP} and modified for use in R-Instat.
#' The modifications include changes to data input and output processes while preserving the core calculation methods.
#' The exponential distribution is commonly used to model the time between events in a Poisson process.
#' The parameters `n`, `r`, and `s` determine the shape and characteristics of the distribution.
#' 
#' @examples
#' # Generate an exponential distribution with minimum value 0, rate parameter 1, and scale 2
#' # x <- seq(0, 5, by = 0.1)
#' # y <- Exponential(x, 0, 1, 2)
#' # plot(x, y, type = "l", main = "Exponential Distribution")
#'
#' @references 
#' For the original source code and more information, please refer to: \href{https://github.com/ET-NCMP/NCMP}{ET-NCMP/NCMP}
#' 
Exponential <- function (x, n, r, s) { (s - n)*(1 - exp(-x / r)) + n }
