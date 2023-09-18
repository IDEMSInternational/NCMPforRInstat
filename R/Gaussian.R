#' @title Gaussian Function
#' 
#' @description
#' This function calculates a Gaussian (normal) distribution with parameters n, r, and s, 
#' for a given input vector x.
#' 
#' @param x A numeric vector or scalar containing input values.
#' @param n The center of the Gaussian curve.
#' @param r The standard deviation (width) of the Gaussian curve.
#' @param s The scale (height) of the Gaussian curve.
#'
#' @return
#' A numeric vector of the same length as x, representing the Gaussian distribution.
#' 
#' @export
#' 
#' @details
#' This function is adapted from the \href{https://github.com/ET-NCMP/NCMP}{ET-NCMP/NCMP} and modified for use in R-Instat.
#' The modifications include changes to data input and output processes while preserving the core calculation methods.
#' The Gaussian function is used to model a normal distribution with parameters n, r, and s. 
#' It describes the probability distribution of a continuous random variable.
#'
#' @examples
#' # Generate a Gaussian curve with center at 0, standard deviation 1, and scale 2
#' # x <- seq(-3, 3, by = 0.1)
#' # y <- Gaussian(x, 0, 1, 2)
#' # plot(x, y, type = "l", main = "Gaussian Distribution")
#'
#' @references 
#' For the original source code and more information, please refer to: \href{https://github.com/ET-NCMP/NCMP}{ET-NCMP/NCMP}
#' 
Gaussian <- function (x, n, r, s) { (s - n) * (1 - exp(-(x ^ 2) / (r ^ 2))) + n }




