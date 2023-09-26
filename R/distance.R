#' @title Calculate Distance from Latitude and Longitude
#' 
#' @description
#' This function calculates the great-circle distance between two points on Earth 
#' given their latitude and longitude coordinates.
#'
#' @param lat1 Latitude of the first point in degrees.
#' @param long1 Longitude of the first point in degrees.
#' @param lat2 Latitude of the second point in degrees.
#' @param long2 Longitude of the second point in degrees.
#'
#' @return
#' The distance in kilometers between the two points.
#' 
#' @export
#' 
#' @details
#' This function is adapted from the \href{https://github.com/ET-NCMP/NCMP}{ET-NCMP/NCMP} and modified for use in R-Instat.
#' The modifications include changes to data input and output processes while carefully maintaining consistent calculations 
#' processes as in the original files.This function calculates the great-circle distance using the Haversine formula. 
#' It takes latitude and longitude values in degrees and returns the distance in kilometers.
#'
#' @examples
#' # Calculate the distance between two points
#' # distance(40.7128, -74.0060, 34.0522, -118.2437)
#' 
#' @references 
#' For the original source code and more information, please refer to: \href{https://github.com/ET-NCMP/NCMP}{ET-NCMP/NCMP}
#' 
distance <- function (lat1, long1, lat2, long2) {
  x1 <- lat1 * pi / 180
  x2 <- lat2 * pi / 180
  y1 <- long1 * pi / 180
  y2 <- long2 * pi / 180
  d1 <- sin(x1) * sin(x2)
  d2 <- cos(x1) * cos(x2) * cos(y1 - y2)
  acos(pmin(pmax(d1 + d2, -1.0), 1.0)) * 6371.009 } # acos is arccos, 6371.009 (km) is WGS84 uniform sphere
