## Functions that are not needed in this package as they are not used.
## They are remnants from long ago that were never removed 


## Function to pull out from the interim data the std. summary stats 
## we need for estimating the rates
extractSummaryMeasures <- function( obsEDIobj, enrollPeriod )
{
  
  enrollTime <- ceiling( obsEDIobj$enrollTime )
  nEnroll <- length( enrollTime )
  
  ## Time during which enrollment happened in interim data 
  enrollWeeks <- min( enrollTime) : max( enrollTime )
  wksPerPeriod <- table( findInterval( enrollWeeks, enrollPeriod$start) )
  
  ## which periods are represented in 'cntsByPeriod'
  w.periods <- as.integer( names( wksPerPeriod ) )
  
  relativeRates <- enrollPeriod$relativeRates
  weighted_nWeeksEnroll <- sum( relativeRates[w.periods] * wksPerPeriod )
  
  
  ## dropout info: number of dropouts, person weeks at risk
  nDropouts <- sum( !is.na( obsEDIobj$dropTime ) )
  personWeeksAtRisk_Dropout <- sum( obsEDIobj$obsTime )
  
  
  ## infection info: number infected, person weeks at risk
  ## For infected ppts, measure time at risk through midpoint of
  ## (lastNegTestTime, infecDxTime).  For others through lastNegTestTime.
  ## 
  nInfected <- sum( !is.na( obsEDIobj$infecDxTime ) )
  
  ## average the values: (lastNegTestTime, infecDxTime)
  aveLastNeg_DxTime <-
    rowMeans( obsEDIobj[, c("lastNegTestTime","infecDxTime")], na.rm=TRUE)
  personWeeksAtRisk_Infection <- sum( aveLastNeg_DxTime )
  
  list( enrollment = list( nEnroll, weighted_nWeeksEnroll ),
        dropout = list( nDropouts, personWeeksAtRisk_Dropout ),
        infection = list( nInfected, personWeeksAtRisk_Infection )
  )
}


enrollmentEstimationFunction <- function( obsEDIobj, enrollPeriod )
{
  enrollTime <- ceiling( obsEDIobj$enrollTime )
  nEnroll <- length( enrollTime )
  
  ## Time during which enrollment happened in interim data 
  enrollWeeks <- min( enrollTime) : max( enrollTime )
  wksPerPeriod <- table( findInterval( enrollWeeks, enrollPeriod$start) )
  
  ## which periods are represented in 'cntsByPeriod'
  w.periods <- as.integer( names( wksPerPeriod ) )
  
  relativeRates <- enrollPeriod$relativeRates
  weighted_nWeeksEnroll <- sum( relativeRates[w.periods] * wksPerPeriod )
  
  ## return estimate
  (nEnroll / weighted_nWeeksEnroll)
}

dropoutEstimationFunction <- function( obsEDIobj, enrollPeriod )
{
  ## dropout info: number of dropouts, person weeks at risk
  nDropouts <- sum( !is.na( obsEDIobj$dropTime ) )
  personWeeksAtRisk_Dropout <- sum( obsEDIobj$obsTime )
  
  ## return estimate
  ( nDropouts / personWeeksAtRisk_Dropout )
}


infectionUpdateFunction <- function( obsEDIobj, enrollPeriod, hyperParams )
{
  ## infection info: number infected, person weeks at risk
  ## For infected ppts, measure time at risk through midpoint of
  ## (lastNegTestTime, infecDxTime).  For others through lastNegTestTime.
  ## 
  nInfected <- sum( !is.na( obsEDIobj$infecDxTime ) )
  
  ## average the values: (lastNegTestTime, infecDxTime)
  aveLastNeg_DxTime <-
    rowMeans( obsEDIobj[, c("lastNegTestTime","infecDxTime")], na.rm=TRUE)
  personWeeksAtRisk_Infection <- sum( aveLastNeg_DxTime )
  
  mapply( sum, hyperParams, list(nInfected, personWeeksAtRisk_Infection),
          SIMPLIFY = FALSE )
}



## Function that produces a posterior object that is constant for some
## parameters and not for others.  
##
## If a function is provided in 'estimationFunction' for a parameter,
## then that parameter will be estimated from the data using that function
## and that value will be taken as fixed for PET prediction

## If a function is provided in 'updateFunction' for a parameter,
## then that parameter will be updated from the data using that function
##
## paramEstimationFunction should be a named list, with names matching parameters
constructPartialPosteriorObject <-
  function( obsEDIobj, enrollPeriod, priorObj, 
            estimationFunction=NULL, updateFunction=NULL )
  {
    
    if ( any( names(estimationFunction) %in% names(updateFunction)) )
      stop("Parameters cannot be listed in both estimationFunction and updateFunction\n")
    
    posteriorObj <- priorObj
    
    
    ## define constant function 'Const'
    Const <- function(n, p) {
      if (length(p)==1) return( rep(p, times=n) )
      stop("Constant value must have length==1\n")
    }
    
    ## if there are parameters to be estimated, then do:
    for ( param in names( estimationFunction ) )
    {
      ## estimate value for parameter
      posteriorObj$params[[ param ]] <- 
        estimationFunction[[ param ]]( obsEDIobj, enrollPeriod )
      
      ## assign new function to be constant
      posteriorObj$Function[[ param ]] <- Const 
    }
    
    ## if there are parameters to be updated, then do:
    for ( param in names( updateFunction ) )
    {
      ## update parameter based on observed data
      posteriorObj$params[[ param ]] <- 
        updateFunction[[ param ]]( obsEDIobj, enrollPeriod, 
                                   priorObj$param[[ param ]] )
    }
    
    ## return posterior object
    posteriorObj
  }            



constructPosteriorObject <- function( obsEDIobj, enrollPeriod, priorObj )
{
  
  ## Since we're using conjugate priors, the form of the posterior
  ## is same as the prior (just different param values)
  posteriorFunc <- priorObj$Function
  
  ## extract summary data from interim observed data object 
  updateList <- extractSummaryMeasures( obsEDIobj, enrollPeriod )
  
  posteriorParams <- priorObj$params
  ## Update the priors with observed information 
  for ( comp in names(posteriorParams) ) 
  {
    posteriorParams[[comp]] <- 
      mapply( sum, posteriorParams[[comp]], updateList[[comp]],
              SIMPLIFY = FALSE )
  }
  
  list( Function = posteriorFunc, params = posteriorParams )
}


## Usage:  
## The function is set up so that the user must define 
##   *either* 'nWeeks' (number of weeks of data to simulate)
##   *or* 'nEnroll' (number of partipants to enroll).
##
## The argument 'maxEnrollment' can be used along with 'nWeeks' if
## desired, to put a cap on enrollment.
##
## The function uses three "get" functions (in 'getEDIfuncs.R')
## as well as "visit schedule" functions (in 'visitSchedule.R'),

## For simulating a certain amount of data from the beginning of a trial
simulateObservedEDIdata <-
  function(rateParams, trtAssignProb, 
           infecRates, protocolVisitFunc,
           enrollPeriod, startWeek=1, 
           nEnroll=NULL, nWeeks=NULL,
           maxEnrollment=NULL)
  {
    
    ## very basic sanity check on args 'nEnroll' and 'nWeeks' ##
    if ( is.null(nEnroll) && is.null(nWeeks) ) 
      stop("You must specify either nEnroll or nWeeks\n") 
    
    if ( !is.null(nEnroll) && !is.null(nWeeks) ) 
      stop("You must specify only *one* of nEnroll or nWeeks\n")
    
    ## ---------------- LET THE GAMES BEGIN ----------------
    
    ## get enrollment count and times
    enr <- getEnrollment(
      rate = rateParams$enrollRate, 
      enrollPeriod = enrollPeriod,
      nEnroll = nEnroll,
      nWeeks = nWeeks, 
      startWeek = startWeek,
      maxEnroll= maxEnrollment )
    
    ## extract enrollment count
    N <- enr$N
    
    ## generate dropout and data for the 'N' participants
    dropout <- getDropout(N, rateParams$dropRate)
    
    ## generate treatment assignments
    nTrt <- length(trtAssignProb)
    trtAssignments <- getTreatmentAssignment(N, prob = trtAssignProb, 
                                             blockSize = 10*nTrt )
    
    ## generate infection data for the 'N' participants
    infect <- getInfection(N, baseRate=rateParams$infecRate, 
                           relRates = infecRates, 
                           trtAssgn = trtAssignments)
    
    
    ## If 'nWeeks' was specified, then we've established a bound on our
    ## follow-up time for participants. So we need to censor the dropout
    ## and infection times using that info.
    if ( !is.null( nWeeks ) )
    {
      ## amount of trial time simulated within this invocation of simulateEDI()
      simTime <- (startWeek - 1 + nWeeks) - enr$times
      
      ## censor dropout and infection times by 'simTime'
      dropout[ is.TRUE(dropout > simTime) ] <- NA
      infect[ is.TRUE(infect > simTime) ] <- NA
      
    } else
    {
      simTime <- NA
    }
    
    
    ## Censor infections that occur after dropout, as they can't be observed
    infect[ is.TRUE( infect > dropout ) ] <- NA
    
    
    ## *** NOTE  ****************************************************
    ##
    ##  We still have records with both dropout and infection times
    ##  recorded, and those all have dropout > infect.
    ##
    ##  We need to keep both times until we calculate the "DX time"
    ##  for the infections, then we can compare dropout to the DX 
    ##  time and keep the earlier of the two.
    ## **************************************************************
    
    
    ## ------------------------------------------------------------ ##
    ## Now, we start to construct the "observed" data, by applying  ##
    ## the study visit Map to the simulated data.  This will allow  ##
    ## us to figure out when: (a) the infections are diagnosed,     ##
    ## (b) the last visit occured (during the simulated period).    ##
    ## ------------------------------------------------------------ ##
    
    ## create observed data object to fill in
    obsEDI <- data.frame( trt = trtAssignments,
                          enrollTime = enr$times,
                          dropTime  = dropout,
                          infecDxTime = NA,
                          lastNegTestTime = NA,
                          obsTime = NA )
    
    ## get infection diagnosis dates for each non-NA infection time
    nonNA.inf <- !is.na( infect )
    
    ## Do the following only if there's at least one infection
    ## (NOTE - to reduce complexity, this section has redundant computation)
    if ( !all( is.na(infect) ) )
    {
      infecDX <- getVisitWeek( infect, protocolVisitFunc, whichVisit = "next")
      
      ## Compute the minimum of the infecDX time, dropout time and simTime
      ## Any events occuring *strictly* after this time are censored
      minTime <- pmin( infecDX, dropout, simTime, na.rm=TRUE)
      
      ## censor infecDX and dropout times
      infecDX[ is.TRUE(infecDX > minTime) ] <- NA
      dropout[ is.TRUE(dropout > minTime) ] <- NA
      
      ## store info in obsEDI
      obsEDI$infecDxTime <- infecDX
      obsEDI$dropTime <- dropout
    }
    
    
    ## Now that we're done making changes to 'dropout' and 'infecDxTime',
    ## we compute the amount of "observation time" ('obsTime') for each
    ## ppt.  This is equal to the 'simTime' for ppt.s without events,
    ## and equal to the event time for participants with events.
    obsEDI$obsTime <- pmin( obsEDI$infecDxTime, 
                            obsEDI$dropTime, simTime, na.rm=TRUE)
    
    
    ## Fill in the time of the last Negative HIV test.  This will be the
    ## last visit prior to 'obsTime'.  Why?  Because 'obsTime' equals
    ## infecDxTime for infecteds, dropTime for dropouts, and simTime for
    ## everyone else, and this is what we want.
    obsEDI$lastNegTestTime <- getVisitWeek( obsEDI$obsTime, 
                                            protocolVisitFunc, "previous")
    
    ## return obsEDI
    obsEDI
  }


summInterimData <- function( obsEDI )
{
  ## enrollment info: *weighted* number of enrollees, weeks enrollment
  enrollTime <- ceiling( obsEDI$enrollTime )
  nWeeksEnroll <- max( enrollTime )
  
  ## dropout info: number of dropouts, person weeks at risk
  nDropouts <- sum( !is.na( obsEDI$dropTime ) )
  personWeeksAtRisk_Dropout <- sum( obsEDI$obsTime )
  
  ## infection info: number infected, person weeks at risk
  ## For infected ppts, measure time at risk through midpoint of
  ## (lastNegTestTime, infecDxTime).  For others through lastNegTestTime.
  ## 
  nInfected <- sum( !is.na( obsEDI$infecDxTime ) )
  
  ## average the values: (lastNegTestTime, infecDxTime)
  aveLastNeg_DxTime <- 
    rowMeans( obsEDI[, c("lastNegTestTime","infecDxTime")], na.rm=TRUE)
  personWeeksAtRisk_Infection <- sum( aveLastNeg_DxTime )
  
  cat("nWeeksEnroll =", nWeeksEnroll,
      "nDropouts =", nDropouts,
      "nInfected =", nInfected, "\n")
}


