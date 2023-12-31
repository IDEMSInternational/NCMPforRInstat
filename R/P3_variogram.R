#' @title Calculate Variograms for Climate Indices
#'
#' @description
#' This function computes variograms for climate indices based on station data. 
#' 
#' @param a2 A data frame containing climate indices data.
#' @param station Name of the station column in 'a2'.
#' @param year Name of the year column in 'a2'.
#' @param month Name of the month column in 'a2'.
#' @param nyb Start year for the variogram period.
#' @param nye End year for the variogram period.
#' @param ne (Optional) Specifies which climate indices to compute variograms for (1 to 8).
#' @param station_df A data frame with station metadata.
#' @param name Name of the station column in 'station_df' (to link with 'a2' using 'station').
#' @param lat Name of the latitude column in 'station_df'.
#' @param lon Name of the longitude column in 'station_df'.
#' @param ncmp_folder Path to the folder containing NCMP data.
#' 
#' @return
#' A data frame containing variogram information for each climate index.
#' 
#' @export
#'
#' @details
#' This function is adapted from the \href{https://github.com/ET-NCMP/NCMP}{ET-NCMP/NCMP} and modified for use in R-Instat.
#' The modifications include changes to data input and output processes while preserving the core calculation methods.
#' The function calculates variograms for climate indices, which are used in geostatistical analysis.
#' Variograms quantify spatial variability and help in predicting values at unobserved locations.
#' The function computes variograms for different climate indices specified by 'ne' and returns
#' information about the variogram models fitted to the data.
#' 
#' @examples
#' # Example usage of p3_variogram:
#' # Assuming 'a2', 'station', 'year', 'month', 'station_df', 'name', 'lat', 
#' #            'lon', and 'ncmp_folder' are defined.
#' # result <- p3_variogram(a2, station, year, month, nyb = 1981, nye = 2010, 
#' #  ne = 1:8,station_df, name, lat, lon, ncmp_folder)
#'
#' @references 
#' For the original source code and more information, please refer to: \href{https://github.com/ET-NCMP/NCMP}{ET-NCMP/NCMP}
#'
#'
p3_variogram <- function(a2, station, year, month, nyb = 1981, nye = 2010, ne = 1:8,station_df, name, lat, lon,ncmp_folder) {
  
  ###################################################################################
  # Set variables with key thresholds for estimating variogram functions            #
  # These are values which have been determined by the ET-NCMP for the purposes of  #
  # generating NCMPs and may differ from other standards or guidlines               #
  # There is currently no minimum threshold for the period to estimate the          #
  # variogram, but testing indicates that short sample periods are problematic      #
  # The ET-NCMP recommendation is to use the climatological period                  #
  # ***** DO NOT CHANGE THE VALUE OF THESE VARIABLES *****                          #
  ###################################################################################
  
  yrlo <- 1900L            # earliest possible year for variograms
  yrhi <- as.POSIXlt(Sys.time())$year + 1899L  # latest possible year == current - 1
  # DP These should be the the same as chosen by user in P2_Indices.R.
  # DP TODO Think of a way to have these set by default without using configuration file.
  # yvlo <- clist.P2$nyb    # recommended start year of variogram period
  # yvhi <- clist.P2$nye    # recommended end year of variogram period
  DmaxT <- 3000           # maximum separation (km) for temperature indices
  DmaxP <- 2000           # maximum separation (km) for precipitation indices
  w <- 20                 # bin width (km)
  
  # Number of stations to process
  # Can exclude stations at the end of the list, but arbitrary selection requires
  # editing of 'P2_Station_List.txt'
  # Still suppressing warning messages - about converting strings to integer
  
  stations_data <- as.character(unique(a2[[station]]))
  stations_metadata <- as.character(unique(station_df[[name]]))
  if (!all(stations_data %in% stations_metadata)) stop("Station information not available for all stations that appear in data.")
  stations <- stations_data
  nstn <- length(stations)
  
  # Generate limits for variogram period
  # Allow the full range of years but effectively recommend the climatological period
  
  # DP Need a way to better link this period to the climatology period chosen in P2.
  if (nye <= nyb) stop("Variogram end period must be after start period.")
  
  # Column names of the possible indices to select i.e. ele[ne]
  ele <- c("TMA", "PrAn", "PrA", "SPI", "TX90p", "TN90p", "TX10p", "TN10p") # NCMP index element
  
  month_name_english <- c("January", "February", "March", "April", "May", "June", "July", 
                          "August", "September", "October", "November", "December")
  
  # Output column names (now standard across CSV files)
  cnames <- c(month_name_english, "Annual")
  
  ###################################################################################
  #    Read the modified station list                                               #
  ###################################################################################
  # With the number of stations already set, extract the subset to use
  
  files <- data.frame(Station = station_df[[name]],
                      Lat = station_df[[lat]],
                      Long = station_df[[lon]], stringsAsFactors = FALSE)
  # Ensure stations are in same order as in the data so that calculations match up.
  files <- files[match(stations, files$Station), ]
  Station <- files$Station
  olats <- files$Lat
  olons <- files$Long
  
  ###################################################################################
  #    Calculate Distances:                                                         #
  # Vectorised calculation - generate lat/long for all pairs of stations            #
  # The use of cosine distance will generate warnings                               #
  # - these can be ignored but need to reset the diagonal (same station) to zero    #
  ###################################################################################
  
  lat1 <- matrix(olats, nrow = nstn, ncol = nstn)
  long1 <- matrix(olons, nrow = nstn, ncol = nstn)
  lat2 <- matrix(olats, nrow = nstn, ncol = nstn, byrow = TRUE)
  long2 <- matrix(olons, nrow = nstn, ncol = nstn, byrow = TRUE)
  D <- distance(lat1, long1, lat2, long2)
  diag(D) <- 0
  D <- round(D, 2)
  
  ncmpn <- c(1L,2L,2L,3L,4L,4L,5L,5L)
  
  # Make year column name consistent with what is expected in the calculation.
  a2$Year <- a2[[year]]
  
  a3_folder <- file.path(ncmp_folder, "A3_Variogram")
  dir.create(a3_folder, recursive = TRUE)
  
  # Loop through all chosen indices
  output_data <- list()
  for (e in ne) {
    # Constrain distances and define bins accordingly
    Dmax1 <- round(max(D), -2)  # Maximum distance between stations (rounded to 100km)
    if (is.element(ncmpn[e], 2:3)) Dmax2 <- DmaxP else Dmax2 <- DmaxT
    Dmax <- min(Dmax1, Dmax2)
    # Currently setting last bin distance as half width past Dmax
    nbin <- as.integer(ceiling(Dmax / w))  # Rounding not actually required
    bin <- (1:nbin) * w                    # Max dist of each bin
    Dl <- bin + w / 2                      # Centre of bin
    
    # Determine which bin each station pair belongs to
    
    L <- pmin(as.integer(floor(D / w)) + 1L, nbin)
    dim(L) <- c(nstn, nstn)
    
    ###################################################################################
    #    Read index data for all stations:                                            #
    ###################################################################################
    # Converting PrAn from percentage back to ratio
    # This also apparently works for PrA (e == 3), although it would be preferable
    # if there was something less arbitrary 
    
    NCMP.stn <- vector("list", nstn)
    for (i in 1:nstn) {
      a2_df <- a2 %>% dplyr::filter(.data[[station]] == stations[i])
      I1 <- tidyr::pivot_wider(a2_df, id_cols = tidyselect::all_of(year), names_from = tidyselect::all_of(month),
                               values_from = tidyselect::all_of(ele[e]))
      I1 <- data.frame(I1)
      if (e == 2L) I1[,2:14] <- I1[,2:14] / 100
      if (e == 3L) I1[,2:14] <- I1[,2:14] / 100 # scale PrA so that matrix inversion is clean
      NCMP.stn[[i]] <- I1
    }
    
    ###################################################################################
    #   Calculate differences in indices and assign into appropriate bin              #
    # Now doing this outside months for simplicity and speed                          #
    # Was storing all station pair differences anyway, so memory is not a factor      #
    # Now utilising data.table, supposed to be faster and more efficient              #
    ###################################################################################
    # Set up variables for using data.table with the station pair loop
    
    yrs <- nyb:nye
    cxnames <- paste(cnames, "x", sep=".")   # column names for first station
    cynames <- paste(cnames, "y", sep=".")   # column names for second station
    Y <- vector("list", nbin)                # empty list for bins
    Idum <- data.table::as.data.table(I1)[0] # dummy data.table
    for (i in 1:nbin) Y[[i]] <- Idum         # empty data.table for each bin
    
    # Loop through pairs of stations
    # For each pair, retain overlapping years only, then cut to variogram period
    # Can then calculate differences and add to the appropriate bin
    for (i in 1:(nstn - 1)) {
      I1 <- data.table::data.table(NCMP.stn[[i]], key="Year")
      for (j in (i + 1):nstn) {
        I12 <- merge(I1, NCMP.stn[[j]])
        I12v <- I12[,.SD[is.element(Year,yrs)]]
        I12vd <- I12v[,cxnames,with = FALSE] - I12v[, cynames, with = FALSE]
        names(I12vd) <- cnames  # consistency with bin data.table
        Y[[L[i,j]]] <- rbind(Y[[L[i,j]]], cbind(I12v[,list(Year)], I12vd))
      }
    } 
    nrec <- sapply(Y, nrow)  # Number of records in each bin 
    
    # Set up Variogram output table for all months
    # This is the information that will be written to file for kriging
    
    X <- data.frame(Month = cnames, Function = NA_character_,
                    n = NA_real_, r = NA_real_, s = NA_real_, "Mean Sq Err" = NA_real_,
                    check.names = FALSE, stringsAsFactors = FALSE)
    Graph <- c("Gaussian", "Exponential", "Spherical")  # Which variogram function to use
    
    ###################################################################################
    # Begins loop over months                                                         #
    ###################################################################################
    # Set specifications for graphs
    # These will be plotted on a single page - have modified the margins
    
    grDevices::pdf(file.path(a3_folder, paste("NCMP", ele[e], "Graph.pdf")))
    graphics::par(mfrow=c(5,3),mar=c(2.5,4,1.5,0.5)+0.1,cex=0.5,mgp=c(1.5,0.4,0),tcl=-0.2)
    
    for (nm in 1:13) {
      
      ###################################################################################
      # Calculate mean of binned values for variogram fitting                           #
      ###################################################################################
      # No need to do this if do not have any data in the bin, but account for missing values
      # Have commented out rounding of the mean binned value
      
      Bl <- rep(NA, nbin)  # empty vector for bin avgs     
      for (k in which(nrec > 0L)) {
        b2 <- Y[[k]][,nm+1,with = FALSE] ^ 2   # Data^2 for this bin and month
        b2 <- b2[[1]]
        Bl[k] <- mean(b2, na.rm = TRUE)       # assign into vector Bl
      }
      #  Bl <- round(Bl,3)
      
      ###################################################################################
      # Fit functional variogram to plot                                                #
      # There are a number of options which could still be considered at this point     #
      # However, testing has found that the major stability issues are                  #
      # too few years (< 20?) or too few (< 10?) stations                               #
      # The current configuration solves for n, r and s from plausible start values     #
      # and apply an iterative loop to set the maximum range <= Dmax                    #
      ###################################################################################
      
      # Set up start values for each parameter
      # s0 == the highest of the mean value of 60-80 and 80-100 percentile bins
      # TEST: Define n0 from binned data - minimum over 0-20 percentile bins
      # TEST: Define s0 and r0 dependant on function to fit
      # TEST: Use 'port' algorithm and constrain parameters
      #     - This is a lot less prone to warnings, and seemingly as plausible
      
      n0 <- max(0.01, 0.95 * min(Bl[1:as.integer(.2 * nbin)], na.rm = TRUE)) 
      sm <- max(Bl[as.integer(.8 * nbin):nbin], na.rm = TRUE)
      s1 <- stats::median(Bl[as.integer(.6 * nbin):(as.integer(.8 * nbin) - 1L)], na.rm = TRUE)
      s2 <- stats::median(Bl[as.integer(.8 * nbin):nbin], na.rm = TRUE)
      s0 <- stats::median(Bl[as.integer(.6 * nbin):nbin], na.rm = TRUE)
      plower <- c(r = 100, s = 0.5 * s0)
      pupper <- c(r = 3 * Dmax, s = sm)
      
      # Gaussian
      clist <- list(r = Dmax, s = sm)
      mod1 <- stats::nls(Bl ~ Gaussian(Dl, 0, r, s), start = clist,
                  algorithm = "port", lower = plower, upper = pupper,
                  control = stats::nls.control(warnOnly = TRUE, maxiter = 100))
      cmod1 <- stats::coef(mod1)
      
      # Exponential - same start values as Gaussian
      mod2 <- stats::nls(Bl ~ Exponential(Dl, 0, r, s), start = clist,
                  algorithm = "port", lower = plower, upper = pupper,
                  control = stats::nls.control(warnOnly = TRUE, maxiter = 100))
      cmod2 <- stats::coef(mod2)
      
      # Spherical
      clist <- list(r = 0.5 * Dmax, s = s0)
      mod3 <- stats::nls(Bl ~ Spherical(Dl, 0, r, s), start = clist,
                  algorithm = "port", lower = plower, upper = pupper,
                  control = stats::nls.control(warnOnly = TRUE, maxiter = 100))
      cmod3 <- stats::coef(mod3)
      
      # Copy fitted coefficients to n,r and s
      # Only constrain n if not using 'port' algorithm,
      # but should then ensure that r and s are also sensible
      
      n <- c(0, 0, 0)
      r <- c(cmod1[1], cmod2[1], cmod3[1])
      s <- c(cmod1[2], cmod2[2], cmod3[2])
      
      # Calculate mean squared error for each fit
      
      E <- rep(Inf, 3)
      E[1] <- sum((Bl -    Gaussian(Dl, n[1], r[1], s[1])) ^ 2, na.rm = TRUE)
      E[2] <- sum((Bl - Exponential(Dl, n[2], r[2], s[2])) ^ 2, na.rm = TRUE)
      E[3] <- sum((Bl -   Spherical(Dl, n[3], r[3], s[3])) ^ 2, na.rm = TRUE)
      
      # Keep name and parameters of best fit
      
      k <- which.min(E)
      X[nm,"Function"] <- Graph[k]
      X[nm,"n"] <- n[k]
      X[nm,"r"] <- r[k]
      X[nm,"s"] <- s[k]
      X[nm,"Mean Sq Err"] <- E[k]
      
      # Plot binned data and best fit variogram
      plot(Dl, Bl, xlim = c(0, 1.05 * Dmax), xaxs = "i", ylim = c(0.0, 1.05 * max(Bl, na.rm = TRUE)), yaxs = "i",
           xlab = "Distance", ylab = "Diff in Index", col = "Blue")
      graphics::title(cnames[nm], line = 0.5)
      
      # Plot all fits in green (commented out but handy for debugging)
      #  curve(get(Graph[1])(x,n[1],r[1],s[1]),col="green",add=TRUE)
      #  curve(get(Graph[2])(x,n[2],r[2],s[2]),col="green",add=TRUE)
      #  curve(get(Graph[3])(x,n[3],r[3],s[3]),col="green",add=TRUE)
      
      # Highlight best fit in red
      graphics::curve(get(Graph[k])(x, n[k], r[k], s[k]), 0, Dmax, col = "red", add = TRUE)
      
      # Default Spherical in black - but this function may not be the best for the data
      # This is not currently used in the decision process, so is commented out full 
      # implementation of a workable default is an open issue
      #  clist <- list(n=0,r=0.5*Dmax,s=s0)
      #  curve(Spherical(x,clist$n,clist$r,clist$s),0,Dmax,col="black",lwd=0.5,add=TRUE)
      
    }
    grDevices::dev.off()  # Close PDF file
    
    ###################################################################################
    # End loop of months                                                              #
    ###################################################################################
    output_data[[ele[e]]] <- X
  }
  X <- dplyr::bind_rows(output_data, .id = "ncmp_index")
  
  ###################################################################################
  #    Write best fit variogram information                                         #
  ###################################################################################
  
  dy <- date()
  attr(X, "dy") <- dy
  attr(X, "nstn") <- nstn
  attr(X, "nyb") <- nyb
  attr(X, "nye") <- nye
  return(X)
}
