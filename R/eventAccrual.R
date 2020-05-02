#' Tabulate event accrual over time since first enrollment
#' 
#' @param trialData either a list of data frames from \code{\link{simTrial}} (i.e., component \code{trialData} from the output list) or a character string specifying a path to an \code{.RData} file outputted by \code{\link{simTrial}}
#' @param atEvents a numeric vector specifying treatment-pooled event counts for which empirical quantiles of time-to-accrual shall be calculated
#' @param atWeeks a numeric vector specifying time points (in weeks) since first enrollment for which empirical quantiles of treatment-pooled event counts shall be calculated
#' @param prob a numeric value in \eqn{(0, 1)} specifying the probability at which the empirical quantiles across the simulated trials are computed (default is 0.5)
#' @param lagTimeMITT a time point (in weeks). Only events with time-to-event greater than \code{lagTimeMITT} are counted in the MITT column (default is 0).
#' @param lagTimePP a time point (in weeks). If specified, only PP events with time-to-event greater than \code{lagTimePP} are counted in the PP column.
#' @param namePP a character string specifying the name of the column in each data frame in \code{trialData} which indicates membership in the PP cohort (default is "\code{pp1}")
tabEventAccrual <- function(trialData, atEvents=NULL, atWeeks=NULL, prob=0.5, lagTimeMITT=0, lagTimePP=NULL, namePP="pp1"){
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
    # generate a list (across trials) of lists (MITT and possibly PP) of vectors of time points (in weeks) at which the event totals in 'atEvents' are reached
    atEventsTimes <- lapply(trialData, function(data, atEvents, lagTimeMITT, lagTimePP, namePP){
      # observed follow-up time
      data$fuTime <- data$exit - data$entry
      
      # indicator of an MITT event, where MITT is defined by 'lagTimeMITT'
      data$eventMITT <- as.numeric(data$event==1 & data$fuTime>lagTimeMITT)
      
      # indicator of an PP event, where PP is defined by 'lagTimePP' and 'namePP'
      if (!is.null(lagTimePP)){ data$eventPP <- as.numeric(data$event==1 & data$fuTime>lagTimePP & data[, namePP]==1) }
      
      # time of the first enrollment
      t0 <- min(data$entry)
      
      # time since first enrollment until the end of follow-up
      data$timeSinceT0 <- data$exit - t0
      
      # order participants by their exit in calendar time
      data <- data[order(data$timeSinceT0), ]
      
      # MITT event totals over time
      data$cumEventMITT <- cumsum(data$eventMITT)
      
      # PP event totals over time
      if (!is.null(lagTimePP)){ data$cumEventPP <- cumsum(data$eventPP) }
      
      # get positions of first matches of 'atEvents'
      # 'idx' may contain NAs if some event totals in 'atEvents' are not observed in 'data';
      # NAs, if any, are grouped together at the end of the 'idx' vector
      idx <- match(atEvents, data$cumEventMITT)
      
      # initialize the output vector
      timesMITT <- rep(NA, length(atEvents))
      
      # fill in slots for event counts in 'atEvents' that are observed in 'data'
      timesMITT[!is.na(idx)] <- data[na.omit(idx), "timeSinceT0"]
      
      # output list
      out <- list(timesMITT=timesMITT)
      
      if (!is.null(lagTimePP)){
        idx <- match(atEvents, data$cumEventPP)
        timesPP <- rep(NA, length(atEvents))
        timesPP[!is.na(idx)] <- data[na.omit(idx), "timeSinceT0"]
        out$timesPP <- timesPP
      }
      
      return(out)
    }, atEvents=atEvents, lagTimeMITT=lagTimeMITT, lagTimePP=lagTimePP, namePP=namePP)
  } else {
    
  }
}