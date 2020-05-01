#' Tabulate event accrual over time since first enrollment
#' 
#' @param trialData either a list of data frames from \code{\link{simTrial}} (i.e., component \code{trialData} from the output list) or a character string specifying a path to an \code{.RData} file outputted by \code{\link{simTrial}}
#' @param atEvents a numeric vector specifying treatment-pooled event counts for which empirical quantiles of time-to-accrual shall be calculated
#' @param atWeeks a numeric vector specifying time points since first enrollment (in weeks) for which empirical quantiles of treatment-pooled event counts shall be calculated
#' @param prob a numeric value in \eqn{(0, 1)} specifying the probability at which the empirical quantiles across the simulated trials are computed (default is 0.5)
#' @param lagTimeMITT a time point (in weeks). Only events with time-to-event greater than \code{lagTimeMITT} are counted in the MITT column (default is 0).
#' @param lagTimePP a time point (in weeks). If specified, only PP events with time-to-event greater than \code{lagTimePP} are counted in the PP column.
tabEventAccrual <- function(trialData, atEvents=NULL, atWeeks=NULL, prob=0.5, lagTimeMITT=0, lagTimePP=NULL){
  if ((is.null(atEvents) & is.null(atWeeks)) | (!is.null(atEvents) & !is.null(atWeeks))){
    stop("Exactly one of the two arguments 'atEvents' and 'atWeeks' must be specified.")
  }
  
  # load the list of data frames if a file is specified
  if (is.character(trialData)){
    if (file.exists(trialData)){
      load(trialData)
      trialData <- trialObj$trialData
      rm(trialObj)
    } else {
      stop("The file specified in 'trialData' does not exist.")
    }
  }
  
  if (!is.null(atEvents)){
    # generate a list (across trials) of lists (MITT and PP) of vectors of time points (in weeks) at which the event totals in 'atEvents' are reached
    atEventsTimes <- lapply(trialData, function(data, lagTimeMITT, lagTimePP){
      # time of the first enrollment
      t0 <- min(data$entry)
      
      # time since first enrollment until the end of follow-up
      data$timeSinceT0 <- data$exit - t0
      
      # order participants by their exit in calendar time
      data <- data[order(data$timeSinceT0), ]
      
      # event totals over time
      data$cumEvents <- cumsum(data$event)
      
      
    }, lagTimeMITT=lagTimeMITT, lagTimePP=lagTimePP)
  } else {
    
  }
}