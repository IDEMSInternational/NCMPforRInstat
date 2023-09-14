#' @title Calculate Percent Days Diagnostics in the Presence of Limited Data
#' 
#' @description
#' This function calculates percent days diagnostics when there is an insufficient number of years available for analysis.
#' 
#' @param ci A climdexInput object containing climate data.
#' @param freq The frequency of the analysis, either "monthly" or "annual".
#' @param type The type of percentile analysis to perform, one of "TN10p", "TN90p", "TX10p", or "TX90p".
#'
#' @return A numeric vector containing the calculated percent days diagnostics.
#' 
#' @export 
#'
#' @examples
#' # Example usage:
#' # data <- your_climdexInput_object
#' # result <- percent.days.op.threshold.mod(data, freq = "monthly", type = "TN10p")
#' 
#' # The 'result' vector will contain the calculated percent days diagnostics.
#' 
#' @details
#' This function is adapted from the ET-NCMP/NCMP project (https://github.com/ET-NCMP/NCMP) and modified for use in R-Instat.
#' It specifically corrects the calculation of percent days diagnostics in cases where there are too few years of data.
#' Previous versions generated warnings and unrealistic negative bias in the estimated values.
#' The modifications include changes to data input and output processes while preserving the core calculation methods.
#' This function directly calculates percent days diagnostics and circumvents most climdex.pcic wrapper functions for ease of use.
#'
#' @references
#' For the original source code and more information, please refer to: https://github.com/ET-NCMP/NCMP
#' 
#' 
percent.days.op.threshold.mod <- function (ci,freq,type) 
  {
  stopifnot(inherits(ci,"climdexInput"))
  freq <- match.arg(freq,c("monthly","annual"))
  type <- match.arg(type,c("TN10p","TN90p","TX10p","TX90p"))
  fo <- switch(type,
               TN10p = c("tmin","q10","<"), TN90p = c("tmin","q90",">"),
               TX10p = c("tmax","q10","<"), TX90p = c("tmax","q90",">"))
  vn <- fo[1]
  qn <- fo[2]
  f <- match.fun(fo[3])
  dat <- f(ci@data[[vn]], ci@quantiles[[vn]]$outbase[[qn]][ci@jdays])
  inset <- ci@dates >= ci@base.range[1] & ci@dates <= ci@base.range[2]
  if (sum(inset) > 0L && length(ci@dates) >= 720L) {
    temp.base <- ci@data[[vn]][inset]
    years.base <- as.POSIXlt(ci@dates[inset])$year + 1900L
    years.base.range <- range(years.base)
    byrs <- years.base.range[2] - years.base.range[1] + 1L
    base.thresholds <- ci@quantiles[[vn]]$inbase[[qn]]
    bdim <- dim(base.thresholds)
    dim(base.thresholds) <- c(bdim[1] * bdim[2], bdim[3])
    yday.byr.indices <- ci@jdays[inset] + (years.base - years.base.range[1]) * bdim[1]
    f.result <- f(rep(temp.base, byrs - 1L), base.thresholds[yday.byr.indices,1:(byrs-1)])
    dim(f.result) <- c(length(yday.byr.indices), byrs - 1L)
    dat[inset] <- rowSums(f.result, na.rm = TRUE)/(byrs - 1)
  }
  dat[is.nan(dat)] <- NA
  na.mask <- ifelse(climdex.pcic:::tapply.fast(dat,ci@date.factors[[freq]],
                                               function(x) sum(is.na(x)) > ci@max.missing.days[freq]),NA,1)
  ret <- climdex.pcic:::tapply.fast(dat, ci@date.factors[[freq]], mean, na.rm = TRUE) * 
    100 * na.mask
  ret[is.nan(ret)] <- NA
  ret * ci@namasks[[freq]][[vn]]
  }

