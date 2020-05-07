# quantile2() is the same as quantile() except it returns NA if the proportion of NAs exceeds 'na.ub'; otherwise, NAs are removed before the quantile is computed
quantile2 <- function(x, probs, na.ub=0.2){
  if (mean(is.na(x))>na.ub){
    return(NA)
  } else {
    return(quantile(x, probs=probs, na.rm=TRUE))
  }
}

#' Tabulate event accrual over time since first enrollment
#' 
#' Tabulates side-by-side different event totals and time periods since first enrollment required to accrue the event totals. The user specifies a vector of either event totals or time points since first enrollment, and \code{tabEventAccrual} completes the table.
#' 
#' @param trialData either a list of data frames from \code{\link{simTrial}} (i.e., component \code{trialData} from the output list) or a character string specifying a path to an \code{.RData} file outputted by \code{\link{simTrial}}
#' @param atEvents a numeric vector specifying treatment-pooled event counts for which empirical quantiles of time-to-accrual shall be calculated
#' @param atWeeks a numeric vector specifying time points (in weeks) since first enrollment for which empirical quantiles of treatment-pooled event counts shall be calculated
#' @param prob a numeric value in \eqn{(0, 1)} specifying the probability at which the empirical quantiles across the simulated trials are computed (default is 0.5)
#' @param lagTimeMITT a time point (in weeks). Only events with time-to-event \eqn{\ge} \code{lagTimeMITT} are counted in the MITT column (default is 0).
#' @param lagTimePP a time point (in weeks). If specified, only PP events with time-to-event \eqn{\ge} \code{lagTimePP} are counted in the PP column.
#' @param namePP a character string specifying the name of the column in each data frame in \code{trialData} which indicates membership in the PP cohort (default is "\code{pp1}")
#' @param na.ub a numeric value specifying an upper limit on the fraction of simulated trials that do not reach a given event count in \code{atEvents} to still compute the empirical quantile of time-to-accrual. If the fraction of such trials exceeds \code{na.ub}, \code{NA} will be produced.
#' 
#' @details All time variables use week as the unit of time.
#' 
#' If the user specifies \code{atEvents}, time periods since first enrollment are computed that are needed to observe \code{atEvents} MITT and \code{atEvents} PP events.
#' 
#' If the user specifies \code{atWeeks}, MITT and PP event totals are computed that are observed by \code{atWeeks} weeks since first enrollment.
#' 
#' The function inputs a large number of simulated trial data, and the computed variables (time periods or event totals) are empirical quantiles at probability \code{prob} of the sample distributions (of time periods or event totals). Medians are computed by default.
#' 
#' @return A data frame (with at least two columns) of event totals and associated time periods since first enrollment required to accrue the event totals in the MITT cohort. If an PP cohort is specified (via \code{lagTimePP}), a third column is added.
#' 
#' @examples
#' simData <- simTrial(N=c(1000, rep(700, 2)), aveVE=seq(0, 0.4, by=0.2), 
#'                     VEmodel="half", vePeriods=c(1, 27, 79), enrollPeriod=78, 
#'                     enrollPartial=13, enrollPartialRelRate=0.5, dropoutRate=0.05, 
#'                     infecRate=0.04, fuTime=156, 
#'                     visitSchedule=c(0, (13/3)*(1:4), seq(13*6/3, 156, by=13*2/3)),
#'                     missVaccProb=c(0,0.05,0.1,0.15), VEcutoffWeek=26, nTrials=5, 
#'                     blockSize=NULL, stage1=78, randomSeed=300)
#' 
#' ## user specifies MITT event totals
#' tabEventAccrual(simData$trialData, atEvents=seq(10, 100, by=10))
#' 
#' ## user specifies MITT and PP event totals
#' tabEventAccrual(simData$trialData, atEvents=seq(10, 100, by=10), lagTimePP=6)
#' 
#' ## user specifies time points since first enrollment
#' tabEventAccrual(simData$trialData, atWeeks=seq(52, 156, by=8), lagTimePP=6)
#' 
#' @seealso \code{\link(simTrial)}
#' 
#' @export
tabEventAccrual <- function(trialData, atEvents=NULL, atWeeks=NULL, prob=0.5, lagTimeMITT=0, lagTimePP=NULL, namePP="pp1", na.ub=0.2){
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
    atEvents <- sort(atEvents)
    
    # generate a list (across trials) of lists (MITT and possibly PP) of vectors of time points (in weeks) at which the event totals in 'atEvents' are reached
    atEventsTimes <- lapply(trialData, function(data, atEvents, lagTimeMITT, lagTimePP, namePP){
      # observed follow-up time
      data$fuTime <- data$exit - data$entry
      
      # indicator of an MITT event, where MITT is defined by 'lagTimeMITT'
      data$eventMITT <- as.numeric(data$event==1 & data$fuTime>=lagTimeMITT)
      
      # indicator of an PP event, where PP is defined by 'lagTimePP' and 'namePP'
      if (!is.null(lagTimePP)){ data$eventPP <- as.numeric(data$event==1 & data$fuTime>=lagTimePP & data[, namePP]==1) }
      
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
      # 'idx' may contain NAs if some event totals in 'atEvents' are not observed in 'data'
      idx <- match(atEvents, data$cumEventMITT)
      
      timesMITT <- rep(NA, length(atEvents))
      
      # fill in slots for event counts in 'atEvents' that are observed in 'data'
      timesMITT[!is.na(idx)] <- data[na.omit(idx), "timeSinceT0"]
      
      # initialize the output list
      out <- list(timesMITT=timesMITT)
      
      if (!is.null(lagTimePP)){
        idx <- match(atEvents, data$cumEventPP)
        timesPP <- rep(NA, length(atEvents))
        timesPP[!is.na(idx)] <- data[na.omit(idx), "timeSinceT0"]
        out$timesPP <- timesPP
      }
      
      return(out)
    }, atEvents=atEvents, lagTimeMITT=lagTimeMITT, lagTimePP=lagTimePP, namePP=namePP)
    
    # generate a vector of empirical quantiles (at probability 'prob') of the time distributions from the first enrollment until each MITT event count in 'atEvents'
    atEventsTimeQtilesMITT <- round(apply(matrix(sapply(atEventsTimes, "[[", "timesMITT"), nrow=length(atEvents)), 1, quantile2, probs=prob, na.ub=na.ub), digits=1)
    
    # initialize the output data frame
    out <- data.frame(atEvents, atEventsTimeQtilesMITT)
    colnames(out) <- c("events", paste0("weeks_qtile", prob, "_MITT"))
    
    if (!is.null(lagTimePP)){
      # generate a vector of empirical quantiles (at probability 'prob') of the time distributions from the first enrollment until each PP event count in 'atEvents'
      out$atEventsTimeQtilesPP <- round(apply(matrix(sapply(atEventsTimes, "[[", "timesPP"), nrow=length(atEvents)), 1, quantile2, probs=prob, na.ub=na.ub), digits=1)
      colnames(out)[3] <- paste0("weeks_qtile", prob, "_PP")
    }
    
    return(out)
    
  } else {
    atWeeks <- sort(atWeeks)
    
    # generate a list (across trials) of lists (MITT and possibly PP) of vectors of time points (in weeks) at which the event totals in 'atEvents' are reached
    atWeeksEvents <- lapply(trialData, function(data, atWeeks, lagTimeMITT, lagTimePP, namePP){
      # observed follow-up time
      data$fuTime <- data$exit - data$entry
      
      # indicator of an MITT event, where MITT is defined by 'lagTimeMITT'
      data$eventMITT <- as.numeric(data$event==1 & data$fuTime>=lagTimeMITT)
      
      # indicator of an PP event, where PP is defined by 'lagTimePP' and 'namePP'
      if (!is.null(lagTimePP)){ data$eventPP <- as.numeric(data$event==1 & data$fuTime>=lagTimePP & data[, namePP]==1) }
      
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
      
      # get positions of the last exit time points since first enrollment before or at 'atWeeks'
      # 'idx' may contain NAs if time points in 'atWeeks' are before the first exit time point since first enrollment in 'data'
      idx <- sapply(atWeeks, function(atWk){ ifelse(data$timeSinceT0[1]<=atWk, max(which(data$timeSinceT0<=atWk)), NA) })
      
      eventsMITT <- rep(NA, length(atWeeks))
      
      # fill in slots for time points in 'atWeeks' that are after the first exit time sice first enrollment
      eventsMITT[!is.na(idx)] <- data[na.omit(idx), "cumEventMITT"]
      
      # initialize the output list
      out <- list(eventsMITT=eventsMITT)
      
      if (!is.null(lagTimePP)){
        eventsPP <- rep(NA, length(atWeeks))
        eventsPP[!is.na(idx)] <- data[na.omit(idx), "cumEventPP"]
        out$eventsPP <- eventsPP
      }
      
      return(out)
    }, atWeeks=atWeeks, lagTimeMITT=lagTimeMITT, lagTimePP=lagTimePP, namePP=namePP)
    
    # generate a vector of empirical quantiles (at probability 'prob') of the event count distributions up to each time point in 'atWeeks'
    atWeeksEventQtilesMITT <- apply(matrix(sapply(atWeeksEvents, "[[", "eventsMITT"), nrow=length(atWeeks)), 1, quantile2, probs=prob, na.ub=na.ub)
    
    # initialize the output data frame
    out <- data.frame(atWeeks, atWeeksEventQtilesMITT)
    colnames(out) <- c("week", paste0("events_qtile", prob, "_MITT"))
    
    if (!is.null(lagTimePP)){
      # generate a vector of empirical quantiles (at probability 'prob') of the time distributions from the first enrollment until each PP event count in 'atEvents'
      out$atWeeksEventQtilesPP <- apply(matrix(sapply(atWeeksEvents, "[[", "eventsPP"), nrow=length(atWeeks)), 1, quantile2, probs=prob, na.ub=na.ub)
      colnames(out)[3] <- paste0("events_qtile", prob, "_PP")
    }
    
    return(out)
  }
}