# Detects anomalies in a time series using Cyclic hybrid ESD (C-H-ESD).
#
# Anomaly Detection Using Cyclic Hybrid ESD Test ----GM
#
# A technique for detecting anomalies in univariate time
# series where the input is a series of observations.
# @name AnomalyDetection
# @param x Time series as a column data frame, list, or vector,
#  where the column consists of the observations.
# @param max_anoms Maximum number of anomalies that C-H-ESD will
# detect as a percentage of the data. The value can be from 0 to 1.
# @param alpha The level of statistical significance with which
# to accept or reject anomalies.
# @param use_decomp If set to \code{'FALSE'} it gives the possibility 
# to detect outliers with the generalized ESD method on the orginal data.  
# By default is set to \code{'TRUE'} and time series decomposition
# is performed before the analysis.
# @param period Defines the number of observations in a single
# period, and used during seasonal decomposition.
# @param e_value Add an additional column to the anoms output
# containing the expected value.
# @param verbose Additionally printing for debugging.

# @details
# @return The returned value is a list with the following components.
# @return \item{anoms}{Data frame containing index, decomposition components, and
# optionally expected values.}
# @return One can save \code{anoms} to a file in the following fashion:
# \code{write.csv(<return list name>[["anoms"]], file=<filename>)}
# @references Rosner, B., (May 1983), "Percentage Points for a
# Generalized ESD Many-Outlier Procedure", Technometrics, 25(2),
# pp. 165-172.
# @export


anomaly_detection = function(x, max_anoms=0.49, alpha=0.01, decomp = "cffilter", period=2, verbose = FALSE, ModeDeviation = NULL){
  
  # Check for supported inputs types
  if(is.null(period)){
    stop("Period must be set to the number of data points in a single period")
  }
  if(is.vector(x) && is.numeric(x)) {
    x <- ts(x, frequency = period)
  } else if(is.ts(x)) {
  } else {
    stop("data must be a time series object or a vector that holds numeric values.")
  }
  
  # Handle NAs
  if (length(rle(is.na(c(NA,x,NA)))$values)>3){
    stop("Data contains non-leading NAs. We suggest replacing NAs with interpolated values (see na.approx in Zoo package).")
  } else {
    x <- na.omit(x)
  }
  
  # Sanity check all input parameterss
  if(max_anoms > .49){
    stop(paste("max_anoms must be less than 50% of the data points (max_anoms =", round(max_anoms*length(x), 0), " data_points =", length(x),")."))
  }
  
  if(!(0.01 <= alpha & alpha <= 0.1)){
    print("Warning: alpha is the statistical significance level, and is usually between 0.01 and 0.1")
  }
  
  ############## -- Main analysis: Perform C-H-ESD -- #################
  # -- Step 1: Decompose data. This will return two more components: trend and cycle   
  if(decomp == "cffilter"){
    # Christiano-Fitzgerald filter
    x_ts <- ts(x, frequency = period)
    # you could add an extra parameter for pu to smooth the trend line
    x_cf <- cffilter(x_ts ,pl=2,pu=200,type='symmetric')
    ToRemove <- unique(c(which(x == 0), 
                         which(is.na(x_cf$trend)),
                         which(is.na(x_cf$cycle))))
    # start creating the object to retrieve anomalies
    data_dec <- data.frame(index = 1:length(x), 
                           Trend = as.vector(x_cf$trend),
                           Cycle = as.vector(x_cf$cycle))
  } else { 
    # run Loess regression to retrieve Trend line
    df<- data.frame(x = 1:length(x), y =x)
    loess50<-loess(y ~ x, df, span=0.55)
    smooth50 <- predict(loess50) 
    # plot(df$x, df$y, pch=19, main='Loess Regression Models')
    # lines(smooth50, x=df$x, col='red')
    ToRemove <- which(x == 0)
    data_dec <- data.frame(index = 1:length(x), 
                           Trend = smooth50,
                           Cycle = (smooth50 - x))
  }
  if(length(ToRemove)>0)
    data_dec <- data_dec[-ToRemove,]
  
  #### Remove data that are too far from the Mode 
  if(!is.null(ModeDeviation)){
    Trend_unique <- unique(trunc(data_dec$Trend))
    ModeTrend <- Trend_unique[which.max(tabulate(match(trunc(data_dec$Trend), Trend_unique)))]
    # here calculate how many standard deviation from the Trend you should remove the values
    DeviationMode <- sd(data_dec$Trend) * ModeDeviation 
    ToRemoveMode <- c(which(data_dec$Trend  > (ModeTrend + DeviationMode)), which(data_dec$Trend  < (ModeTrend - DeviationMode)))
    if(length(ToRemoveMode)>0){
      data_dec <- data_dec[-ToRemoveMode,]
      ToRemove <- unique(c(ToRemove,ToRemoveMode))
    }
  }
  
  sign_n <- sign( data_dec$Cycle )
  ## to make it compatible to changes in the trend component, I subtract the trend from the original values
  data_dec$Values_Minus_Trend <- (abs( x[data_dec$index] - data_dec$Trend) + abs(data_dec$Cycle)) * sign_n 
  
  
  n <- nrow(data_dec)
  max_outliers <- trunc(n*max_anoms)
  func_mid <- match.fun(median)
  func_dev <- match.fun(mad)
  R_idx <- 1L:max_outliers
  num_anoms <- 0L
  
  
  # Compute test statistic until r=max_outliers values have been 
  # removed from the sample.
  for (i in 1L:max_outliers){
    if(verbose) message(paste(i,"/", max_outliers,"completed"))
    
    # Compute statistics
    CenteredValues = data_dec$Values_Minus_Trend - func_mid(data_dec$Values_Minus_Trend)
    
    Deviation <- func_dev( data_dec$Cycle )
    if(Deviation == 0) 
      break
    
    Scaled <- abs(CenteredValues/Deviation)
    
    R <- max(Scaled)
    
    ## Compute critical value.
    p <- 1 - alpha/(2*(n-i+1))
    
    # Calculate lambda
    t <- qt(p,(n-i-1L))
    lam <- t*(n-i) / sqrt((n-i-1+t**2)*(n-i+1))
    
    temp_max_idx <- which(Scaled == R)[1L]
    R_idx[i] <- data_dec[[1L]][temp_max_idx]
    data_dec <- data_dec[-which(data_dec[[1L]] == R_idx[i]), ]
    
    if(R > lam)
      num_anoms <- i
  }
  
  if(num_anoms > 0) {
    R_idx <- R_idx[1L:num_anoms]
    all_data <- data.frame(index = 1:length(x), anoms = x)
    TotalAnomaliesIndexes <- c(R_idx, ToRemove)
    anoms_data <- subset(all_data, (all_data[[1]] %in% TotalAnomaliesIndexes))
  } else {
    anoms_data <- NULL
  }
  return (list(anoms = anoms_data, num_obs = length(x)))
}

