#' @import stats
NULL
#' @import utils
NULL
#' @import survival
NULL

globalVariables(c("trialObj","futime","event","exit","trt","totInfec","entry","N",
                  "trialListCensor","out","n","Nvacc","approxOk","V","P","pwList",
                  "N.vax.arms", "FR", "FR_loCI", "FR_upCI", "HR", "HR_loCI", "HR_upCI",
                  "VE","VE_loCI","VE_upCI", "evalTime", "nEvents", "nEvents.1", 
                  "nEvents.2", "null.p", "varlogFR"))

## Create a auxillary function "is.TRUE" to use for logical comparisons
## is.TRUE(x) returns TRUE for each element of 'x' that is TRUE, and
##   returns FALSE for all else (i.e. for FALSE and NA)
is.TRUE <- function(x) 
{
  if ( !is.logical(x) ) stop("Argument to 'is.TRUE' must be of type Logical")
  x & !is.na(x)
}


visitScheduleTruncated <-  function(visitWeeks, thruWeek=NULL, nVisits=NULL)
{
  if ( is.null( c(thruWeek, nVisits) ) )
    stop("Must specify one of the arguments: thruWeek or nVisits\n")
  
  ## determine how many more visit weeks need to be added to 'visitWeeks'
  if ( !is.null(nVisits) )
  {
    visitWeeks[ 1:nVisits ]
  } else
  {
    ## find values of 'visitWeeks' that are strictly greater than 'thruWeek'
    gt_thruWeek <- which( visitWeeks > thruWeek )
    
    if ( length(gt_thruWeek) > 0 ) {
      return( visitWeeks[ 1:gt_thruWeek[1] ] )
    } else {
      ## if we get here it means that thruWeek is the last week of our visit schedule
      ## and for some reason we still need to add another week - don't know why...
      ## I'm going to see if my code will work okay if I just return the entire visit
      ## vector and don't add another 'imaginary' visit after the last one. 
      return( visitWeeks )
    }
  }
}

## This function takes as input a vector of times (for 505 we are
## using time in weeks since first vaccination), and returns the
## time of the next or previous scheduled visit.  The returned times 
## are according to the visit schedule, not based off ppt data.
##
## The required argument 'protocolVstFunc', must be a function whose
## first argument takes a time (or vector of times) and that returns
## a vector of all scheduled visit times thru the largest time given
## (e.g. see function 'visitSchedule505')

getVisitWeek <- function( week, visitWeeks, whichVisit=c("next","previous"))
{
  
  whichVisit <- match.arg( whichVisit )
  
  ## check for NAs
  noNAs <- !is.na( week )
  
  wk <- week[ noNAs ]
  
  ## get scheduled visit weeks through max week in 'week'
  schedVstWeeks <- visitScheduleTruncated(visitWeeks, thruWeek=max(wk))
  
  ## if schedVstWeeks has only one value, then need to append
  ## on a 2nd, so avoid errors in 'cut'
  if (length(schedVstWeeks) == 1)
    schedVstWeeks <- c(schedVstWeeks, Inf)
  
  ## find the intervals that the values of 'wk' lie in
  interval <- cut( wk, breaks=schedVstWeeks, right=TRUE, labels=FALSE)
  
  ## If we want the "previous" visit, we return the lower bound of the
  ## interval (schedVstWeeks[interval]), if we want the "next" we return
  ## the upper bound of the interval (schedVstWeeks[interval+1]  )
  if (whichVisit == "previous")
  {
    week[ noNAs ] <- schedVstWeeks[ interval ]
  } else if (whichVisit == "next")
  {
    week[ noNAs ] <- schedVstWeeks[ interval + 1 ]
  }
  
  week 
}

## Function returning a list with components: enrollRate, dropRate, infecRate
## sampled from their joint prior distribution. 
##
## Right now we're assuming independence so we sample from each of three
## prior separately

## This function has a generic interface that must be maintained.
## However the guts of it can be altered as desired as can be the form
## that is assumed for 'paramList'
##
## n = number of sets of parameter values to sample
GammaDist <- function(n, paramList)
{
  do.call( rgamma, as.list( c(n=n, paramList) ) )
}


## function for inputing ...
gammaInput <- function(alpha, beta)
{ 
  list( shape = alpha, rate = beta )
}

## This prior distribution returns the same values passed to it (i.e., it is a non-random point-mass distribution)
Constant <- function(n, paramList)
{
  rep( paramList, times=n )
}

## 'sampleRates' samples 'n' values from each of the probability distributions specified in the 'from' list
## 'from' contains two sublists: 'Function' (distributions) and 'params' (parameters of the distributions)
sampleRates <- function(n, from)
{
  Function <- from$Function # a list
  Params <- from$params     # a list
  
  out <- lapply( 1:length(Params), function(i, n) Function[[i]](n, Params[[i]]), n=n ) 
  names( out ) <- names( Params )
  
  ## rearrange (if needed)
  out <- out[ c("enrollment", "dropout", "infection") ]
  
  ## now rename
  names( out ) <- c("enrollRate", "dropRate", "infecRate")
  
  out
}

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

## Generate enrollment times
## --------------------------
## Sub-function used in 'getEnrollment'.  This specifies the 
## distribution from which the enrollment data come.  This is
## where the actual sampling occurs, and the function to swap out
## if you want to make different distributional assumptions.

generateEnrollment <- function(start, end, rate)
{
  ## number of enrollees during period is poisson distributed with 
  ## rate = 'rate' * (end - start + 1), and the times are uniformly
  ## distributed on the interval (start,end)
  N <- rpois(1, lambda = rate * (end - (start-1)))
  list( N = N,  times = sort( runif(N, min=start-1, max=end) ) )
}


## Generate dropout times
## ----------------------
## Specify the number of participants to generate times for (N)
## and the dropout rate (rate) 
##
## Returns continuously valued times

generateDropout <- function(N, rate)
{
  ## changing the 'rate' parameter to a form that gives the
  ## prob(event in interval (0,1)) = rate, i.e., P(dropout time <=1) = 1 - exp(-lambda) = rate
  rexp(N, rate= -log( 1 - rate) )
}


## Generate infection times  
## ------------------------
## These should be taken as the time at which the infection is
## diagnosable (by Elisa), not when infection occurs.
## 
## Specify the number of participants to generate times for (N)
## and the infection rate (rate) 
##
## Returns continuously valued times

generateInfection <- function(N, rate)
{
  ## the argument 'rate' must either be a numeric vector of length 1 
  ## or a data.frame containing columns: start, end and rate
  if ( length(rate) == 1 && !is.data.frame(rate)) {
    
    ## changing the 'rate' parameter to a form that gives the
    ## prob(event in interval (0,1)) = rate
    rexp(N, rate= -log( 1 - rate) )
  } else
  {
    if ( is.data.frame(rate) )    
    {
      ## ensure the rows of data frame are ordered properly w.r.t time
      rate <- rate[ order( rate$end ), ]
      
      ## for each row of 'rate' generate 'N' infection times.  The final
      ## value will be the minimum of these times (or NA) - see below
      ## for details
      nRows <- nrow(rate)
      offsetTimes <- c(rate$start[1]-1, rate$end[ -nRows ]) 
      rowTimes <- 
        lapply(1:nRows, function(i, N, r, offset, max) {
          
          ## generate exponential with approp. rate and add offset
          times <- offset[i] + rexp(N, -log( 1 - r[i]) )
          
          ## values greater than 'max' get set to missing
          times[ times > max[i] ] <- NA 
          
          times 
        }, N=N, r=rate$rate, offset=offsetTimes, max=rate$end)
      
      do.call(pmin, c(rowTimes, na.rm=TRUE))
      
    } else
      stop("Argument 'rate' must be a data.frame containing both rate ",
           "and time period information in order to utilize multiple ",
           "rates\n")
  }
}

## we typically supply 'rate', 'enrollPeriod' data.frame and 'nEnroll', and keep 'nWeeks' and 'maxEnroll' =NULL

getEnrollment <-
  function( rate, enrollPeriod=NULL, nEnroll=NULL, nWeeks=NULL, 
            startWeek=1, maxEnroll=NULL)
  {
    if ( is.null( enrollPeriod ) )
      enrollPeriod <- data.frame( start=1, end=NA, relativeRate=1 )
    
    nPeriods <- nrow( enrollPeriod )
    
    ## all periods have defined 'start's, but not necessarily 'end's
    startPeriod <- max( which( startWeek >=  enrollPeriod$start ) )
    
    ## define 'endWeek'
    endWeek <- ifelse( !is.null(nWeeks), startWeek + nWeeks - 1, NA)
    
    ## use 'endWeek' if it's not NA, else use max period
    endPeriod <- ifelse( !is.na(endWeek), 
                         max( which(endWeek >= enrollPeriod$start) ),
                         nPeriods)
    
    ## if a single rate was specified, repeat it so length = nPeriods
    if ( length(rate) == 1 )
      rate <- rep(rate, nPeriods)
    
    ## initialize 'N' (number enrolled) and 'times' (enrollment times)
    N <- 0
    times <- NULL
    
    ## if 'nWeeks' is specified, generate all nWeeks of data,
    ## then truncate to 'maxEnroll' if necessary
    if ( !is.null(nWeeks) )
    {
      for ( i in startPeriod:endPeriod )
      {
        ## get weekly enrollment rate for period 'i'
        rate.i <- rate[i] * enrollPeriod$relativeRate[i] 
        
        start.i <- ifelse(i > startPeriod, enrollPeriod$start[i], startWeek)
        end.i <- ifelse(i < endPeriod, enrollPeriod$end[i], endWeek)
        
        ## generateEnrollment returns a list with components 'N' and 'times'
        out <- generateEnrollment(start.i, end.i, rate.i)
        
        N <- sum(N, out$N)
        times <- c( times, out$times )
      }
      if ( !is.null(maxEnroll) && N > maxEnroll )
      {
        N <- maxEnroll
        times <- times[ 1:N ]
      }
    } else {
      ## If 'nEnroll' is specified, then generate data until goal is reached
      for ( i in startPeriod:endPeriod )
      {
        ## get weekly enrollment rate for period 'i'
        rate.i <- rate[i] * enrollPeriod$relativeRate[i] 
        
        start.i <- ifelse(i > startPeriod, enrollPeriod$start[i], startWeek)
        end.i <- ifelse(i < endPeriod, enrollPeriod$end[i], endWeek)
        
        ## If end.i is 'NA' (i.e., not specified), then we set it to a value
        ## much larger than should be needed to generate remaining enrollees
        if ( is.na(end.i) )
          end.i <- start.i + ceiling( 10 * (nEnroll - N)/rate.i )
        
        ## generateEnrollment returns a list with components 'N' and 'times'
        out <- generateEnrollment(start.i, end.i, rate.i)
        
        N <- sum(N, out$N)
        times <- c( times, out$times )
        
        if ( N >= nEnroll ) break 
        
        if ( i == endPeriod )
          stop("Probable error in generateEnrollment(): ",
               " value of 'nEnroll' not reached\n\n")
      } 
      
      ## restrict to the first nEnroll enrollees
      N <- nEnroll
      times <- times[1:N]
    }
    
    list( N = N, times = times)
  }

getDropout <- function(N, rate)
{
  generateDropout(N, rate)
}

getInfection <- function(N, baseRate, relRates=NULL, trtAssgn=NULL)
{
  ## if no object is provided for relRates, then 'trtAssgn' isn't
  ## needed and all data is generated using the 'baseRate' value
  if ( is.null(relRates) )
  {
    ## need to fix this
    return( generateInfection(N, baseRate) )
    
  } else {
    
    nTrts  <- length( unique(relRates$trt) )
    nRates <- length( unique(relRates$relRate) )
    
    ## compute absolute rates from base-rate and relative rates
    relRates$rate <- relRates$relRate * baseRate
    
    ## If only one Trt or only one rate, then don't need treatment info
    if ( nTrts == 1 || nRates == 1 ) {
      ## get infection times for all ppt.s at once
      n <- N
      if (nRates == 1) {
        ## just one rate
        return( generateInfection(N, unique(relRates$rate) ) )
      } else {
        ## just one treatment group (and multiple rates)
        return( generateInfection(N, relRates[,c("rate","start","end")]) )
      }
    } else {
      
      ## More than one treatment and more than one rate - need treatment
      ## assignment info for this case
      if ( is.null(trtAssgn) )
        stop("Treatment assignment information is needed for simulation",
             "of infection times.\n", "Please provide this information ",
             "via argument 'trtAssgn'\n")
      
      ## create numeric vector to store infection times in, since we aren't
      ## able to generate them all at once and simply return them 
      infecTimes <- numeric(N) 
      
      ## get counts for all the treatments
      trtCnt <- table(trtAssgn)
      
      for ( trt in names(trtCnt) ) {
        ## number of ppt.s with treatment 'trt'
        n.trt <- trtCnt[ trt ]
        
        ## identify which ppt.s have treatment 'trt'
        w.trt <- which( trtAssgn == trt )
        
        ## pull out rates just for 'trt'
        rates.trt <- relRates[ relRates$trt == trt, c("rate","start","end")]
        
        ## generate infection times for treatment 'trt' and store
        infecTimes[ w.trt ] <- generateInfection(n.trt, rates.trt)
      }                
      return( infecTimes )
    }    
  }
}


## Function to compute the Greatest Common Divisor of a set of integers
## (used by function getBlockSize() )
##
##  Arguments:  
##      nvec:  a vector of integers
##             The argument value need not be of type integer - it will be coerced
##             to integer internally
##             
##  Return Value:
##      the greatest common divisor (integer)
##
gcd <- function(nvec) {

    nComp <- length(nvec)
    ## Ensure that nvec has at least 2 components
    if ( nComp < 2 )
        stop("Error(gcd): Argument 'nvec' must have length >= 2\n")

    ## ensure that the input vector contains only integers
    if ( any( abs(round(nvec) - nvec) > sqrt(.Machine$double.eps)) )
        stop("Error(gcd): Non-integer values found in input vector\n")

    ## coerce to integer (to avoid dealing with numerical fuzz in comparisons
    nvec <- as.integer( round(nvec) )

    ## Defines a recursive function that computes the gcd of a pair of integers
    ## (taken off an R-mailing list post).
    gcd2 <- function(a,b) ifelse(b==0, a, gcd2(b, a %% b) )

    ## get pairwise gcd between current gcd and the next component
    for (i in 1:(nComp-1)) {
        nvec[i+1] <- gcd2(nvec[i], nvec[i+1])
    }
    nvec[ nComp ]
}


## getBlockSize:  
## A function that computes possible block sizes to use in a blocked randomization.
## It either reports the minimum possible block size (any multiple of it also is
## a usable value), or the smallest block size that falls within a given range
## (e.g. [10, 30] ) if the user provides a value via the 'range' argument.
## If 'range' is specified and there is no block size within that range, the function
## returns NA.
##
##  Arguments:  
##      nvec:  a vector of counts of ppt.s allocated to various treatment arms
##     range:  a length 2 integer vector of form [min, max] that gives the 
##             upper and lower bounds of a region of values to consider for the
##             blockSize of a randomization (the bounds are considered part of the
##             region).  If 'range' is specified then the function returns the
##             smallest multiple of the minimum block size that lies in the region
##             (if any) - or NA if none do.
##             

#' Determine block size for use in blocked randomization
#'
#' \code{getBlockSize} returns the minimum block size (possibly within a specified range) that is compatible with a trial's overall treatment assignment totals.
#' 
#' @param nvec vector specifying the number of participants to be assigned to each treatment group.  The vector should have one component per group, so that its length equals number of groups. The sum of \code{nvec} should equal the total enrollment for the trial.
#' @param range (Optional) vector of length two giving the lower and upper bounds (respectively) on block sizes that the user wishes to consider.
#' 
#' @details The ordering of the components of \code{nvec} is not important, so using \code{nvec = c(x,y,z)} will produce the same results as using \code{nvec = c(z,x,y)}.
#' 
#' In block randomization one does not necessarily want the smallest block size, which is the reason for the existance of the \code{range} argument.  For example, a trial with a 1:1 randomization allocation between two groups would have a minimum block size of 2, which most people would consider to be too small.  So a typical usage of \code{getBlockSize} would be to use \code{range} to set a minimum acceptable block size, through use of vector of form \code{c(lowerBound, Inf)}.  A large trial should probably have a block size on the order of 10-20 or larger, depending on factors including the total trial size and speed of enrollment, so setting a minimum is a good idea.
#'
#' @return An integer or NA.  If the user does not specify \code{range}, then the function will always return an integer, which is the smallest block size compatible with the specified vector of treatment group sizes.  If the user \emph{has} specified the \code{range}, then the function adds the further constraint that the block size must lie in the closed interval given by \code{range} (i.e., the block size must be greater-than-or-equal-to \code{range[1]} and less-than-or-equal-to \code{range[2]}).  If there are no compatible block sizes that lie in the given interval, then an NA is returned.
#' 
#' Note that the value returned is the \strong{minimum} block size that is compatible, not necessarily the only one. Any other compatible block sizes (if any exist) will be integer multiples of the minimum size.  You can check the feasibility of various integer multiples by seeing if they divide evenly into the total trial size (i.e., into the sum of \code{nvec}).
#'
#' @examples 
#'
#' getBlockSize(nvec = c(375, 375) ) 
#' ## specify a minimum block size of 10 (no maximum)
#' getBlockSize(nvec = c(375, 375), range = c(10, Inf) ) 
#' 
#' getBlockSize( nvec = c(30, 510, 390) )
#' ## require a minimum block size of 10 and maximum of 30 
#' ## (not possible with this nvec, so function returns NA)
#' getBlockSize( nvec = c(30, 510, 390), range = c(10, 30) )
#' 
#' @export
getBlockSize <- function(nvec, range=c(0,Inf)) {

    minBlk <- sum( round(nvec) / gcd(nvec) )

    ## if 'range' was specified by user do some checking on the value
    if ( !missing(range) ) {
      if (length(range) != 2)
          stop("Error(getBlockSize): The 'range' parameter must be a vector",
               " of length 2\n")

      if (is.infinite(range[1]))
          stop("Error(getBlockSize): The first value of the 'range' vector is",
               " infinite - this is not allowed\n")

      if ( range[2] < range[1] )
          stop("Error(getBlockSize): The first value of the 'range' vector must",
               " be less-than-or-equal-to the second value\n")

      ## If the range minimum was set to less than 1, reset to 1.  The blockSize
      ## must be an integer and a value below 1 is nonsensical
      if (range[1] < 1 )  range[1] <- 1

      ## Minimum multiple of minBlk needed to meet/exceed the minimum range value
      minMult <- ceiling(range[1]/minBlk)
      ## Maximum multiple of minBlk usable to meet/fall-below the maximum range value
      maxMult <-   floor(range[2]/minBlk)

      ## if the minimum is larger than the maximum, then no value lies in the range
      if ( minMult > maxMult )
        return (NA)
      else
        return ( minMult * minBlk )
    } else {
        ## if range was not specified
        return (minBlk)
    }
}


## Treatment assignment ties in with infection, since infections rates
## will vary across treatments - unless all are equally (in)effective.
getTreatmentAssignment <- 
  function(n, prob=NULL, nPerTrt=NULL, blockSize=NULL, seed=NULL, approxOK=FALSE)
  {
    ## Arguments:  
    ##-----------
    ##   n     (integer) total number of assignments to generate
    ##
    ##   prob  (numeric vector) A (possibly named) vector of assignment
    ##         probabilities for all treatments.  If the elements are named,
    ##         the names are used in the returned object (see details in
    ##         "Return Object" below)
    ##
    ##   nPerTrt  Controls how many ppts are assigned to each trt group.
    ##            Valid values are:
    ##              (character) "fixed" - 
    ##                  Assigns round(n * prob[i]) ppts to the i-th treatment
    ##                  (adjusted if sum( round(n * prob) ) != n ). 
    ##              (character) "random" -
    ##                  Number assigned to each trt is random - distributed as
    ##                  multinomial with probabilities given by 'prob'
    ##              (integer vector)
    ##                  Direct specification of the number ppt.s for each trt
    ##
    ##   blockSize  (integer) Specifies that a blocked randomization approach 
    ##              to treatment assignment should be used, with block size
    ##              size as given by this argument. 
    ##
    ##   seed  (integer) <Optional> Random number seed that can be specified 
    ##         in order to produce reproducible results 
    ##
    ##   approxOK (logical) when nPerTrt=="fixed" is it okay for (prob * n) to
    ##                      not be integers (i.e. is it okay for the number of
    ##                      of ppts in each trt to approximately match 'prob'
    ##                      (rather than exactly)?
    ##
    ## Return Object:
    ## --------------
    ##   A factor vector of length 'n' containing trt assignments. 
    ##   If the argument 'prob' is a named vector then those names are used as
    ##   the levels of the factor that is returned.  Otherwise some default names
    ##   are assigned by attaching a prefix to the treatment's integer code.
    ##   Treatment codes are 1, 2, ..., length(prob), with integer i representing
    ##   treatment with assignment probability of prob[i]
    
    ##
    if ( !is.null(blockSize) ) 
      useBlock <- TRUE
    
    if ( !is.null(prob) ) {
      useProb <- TRUE
      if ( !is.numeric(prob) || any(prob<0 | prob>1) || 
             abs( sum(prob) - 1 ) > .Machine$double.eps ^ 0.5 )
        stop("Argument 'prob' must be a vector of probabilites summing to 1.\n")
    }
    
    if ( !is.null(nPerTrt) ) {
      useNPT <- TRUE
      if ( is.character(nPerTrt) && is.null( prob ) ) 
        stop("Argument 'prob' must be specified too, when 'nPerTrt' is ",
             "set to one of: ('fixed', 'random')\n") 
    }
    
    ## set seed if user provided one
    if ( !is.null(seed) && is.numeric(seed) ) 
      set.seed(seed)
    
    ## if block randomization chosen, do this:
    if ( useBlock ) {
      
      ## the argument 'prob' takes precedence if both it and nPerTrt were given
      if ( useProb ) {  

        ## number of each trt group in each block
        nPerBlock <- round( blockSize * prob )

        ## Verify that the 'blockSize' specified is compatible with the treatment
        ## prob.s given by 'prob'
        if ( sum( nPerBlock ) != blockSize ) {
          stop("\n  The given value of blockSize, ", blockSize,
               ", is not compatible with the specified treatment\n",
               "distribution.  Please use function 'getBlockSize' ",
               "to determine a suitable\n block size (see ",
               "help(getBlockSize) ), and then re-run\n\n")
        }

        ## We allow user to use a poor block size if they really want to (may be
        ## necessary in odd cases?), but issue a warning about it
        if ( max( abs( nPerBlock/sum(nPerBlock) - prob ) > sqrt(.Machine$double.eps) )) {
           warning("\n",
             "The given value of blockSize, ", blockSize, "is not compatible",
             "with the value\n", "of argument 'prob'.  Blocks of the specified ",
             "size cannot create a\n", "treatment distribution that matches ",
             "'probs'.  We suggest using\n", "function 'getBlockSize' to ",
             "determine a suitable block size and then\n", "re-running",
             " (see help(getBlockSize))\n\n", immediate=TRUE)
        } 
      } else {
        if ( useNPT ) {  
          if ( sum(useNPT) != blockSize )
            stop("The values of 'nPerTrt' must sum to the given blockSize \n")
          nPerBlock <- nPerTrt
        } else { 
          stop("Must specify either 'prob' or 'nPerTrt'\n")
        }
      }
      
      ## do Block randomization and return
      nBlocks <- ceiling( n / blockSize)
      
      nTrt <- length(nPerBlock)
      
      ## vector of length 'blockSize' containing the correct number of each trt in it (Doug's version: 'times=nPerBlock')
      trtBlock <- rep(1:nTrt, times=nPerBlock)
      
      ## generates a vector of randomized blocks
      trtVec <- unlist( lapply(1:nBlocks, function(i, trtBlk, bS)
        trtBlk[ order( runif(bS) ) ],
                               trtBlk=trtBlock, bS=blockSize) )
      
      ## keep the first 'n' of them (in case there are more than 'n')
      if ( n < length(trtVec) )
        trtVec <- trtVec[1:n]
      
    } else {   
      
      ## block randomization NOT chosen:
      if ( useProb ) {
        ## 'prob' given and 'nPerTrt' one of "fixed" or "random"
        if ( !useNPT || !nPerTrt %in% c("fixed","random") ) 
          stop("Either 'fixed' or 'random' must be specified for 'nPerTrt' ",
               "when 'prob' is given and blockSize is not.\n") 
        
        nTrt <- length(prob)
        
        if ( nPerTrt == "random" ) {
          trtVec <- sample.int(n = nTrt, size = n, prob = prob, replace=TRUE)
        } else { ## nPerTrt == "fixed"
          np <- round( prob * n ) 
          if ( !approxOk && !all.equal( np, prob * n) )
            stop("The given values of 'n' and 'prob' don't produce an integer",
                 "number of ppt.s for each treatment\n")
          if ( sum(np) != n )
            stop("The sum of round('n'*'prob') must equal'n' \n")
          
          ## create a vector with the correct number of each trt in it, and
          ## then randomize it
          trtVec <- rep.int(1:nTrt, times=np)[ order(runif(n)) ]
        }
      } else {
        ## 'prob' not given, so 'nPerTrt' must be an integer vector
        if ( !all.equal( nPerTrt, as.integer(nPerTrt)) || sum(nPerTrt) != n )
          stop("Value of 'nPerTrt' is not integer and/or does not sum to 'n' \n") 
        nTrt <- length(prob)
        trtVec <- rep.int(1:nTrt, times=nPerTrt)[ order(runif(n)) ]
      }
    }
    
    ## create names for 'prob' if they weren't supplied by user
    if ( !useProb || is.null( names(prob) ) ) {
      prefix <- "Trt"
      trtNames <- paste(prefix, 1:nTrt, sep="")
    } else {
      trtNames <- names(prob)
    }
    
    ## Last step, convert to factor
    trtVec <- factor(trtVec, labels=trtNames)
    
    trtVec
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


## Usage:  
##
## For enrollment, the user must define either:
##   'nWeeksEnroll' (number of weeks of enrollment data to simulate)
##   *or* 'nEnroll' (number of partipants for whom to simulate 
##                   enrollment data).
##
## The argument 'maxEnrollment' can be used along with 'nWeeksEnroll'
## (if desired), to put a cap on enrollment.
##
## The two arguments 'nWeeksFU' and 'nWeeksTrialTime' each specify the 
## amount of follow-up to be performed: 
##   'nWeeksFU' specifies the amount of follow-up for each ppt  
##   'nWeeksTrialTime' specifies total trial duration - each ppt's
##                     follow-up time will vary based on enrollment time
## Only one of these two arguments can be specified.
##
## The function uses three "get" functions (in 'getEDIfuncs.R')
## as well as "visit schedule" functions (in 'visitSchedule.R'),



## For simulating a certain amount of data from the beginning of a trial
simFullEDIdata <-
  function(rateParams, trtAssignProb, blockSize, infecRates, protocolVisits,
           enrollPeriod, startWeek=1, 
           nEnroll=NULL, nWeeksEnroll=NULL, maxEnrollment=NULL,
           nWeeksFU=NULL, nWeeksTrialTime=NULL)
  {
    ## very basic sanity check on arguments 'nEnroll' and 'nWeeksEnroll'
    if ( is.null(nEnroll) && is.null(nWeeksEnroll) ) 
      stop("You must specify either nEnroll or nWeeksEnroll\n") 
    
    if ( !is.null(nEnroll) && !is.null(nWeeksEnroll) ) 
      stop("You must specify only *one* of nEnroll or nWeeksEnroll\n")
    
    if ( !is.null(nWeeksFU) && !is.null(nWeeksTrialTime) ) 
      stop("You must specify only *one* of nWeeksFU or nWeeksTrialTime\n")
    
    if (!is.null(nWeeksEnroll) && !is.null(nWeeksTrialTime) && nWeeksEnroll > nWeeksTrialTime) 
      stop("You have requested ", nWeeksEnroll," weeks of enrollment (via argument",
           " nWeeksEnroll), but\n", "allowed only ", nWeeksTrialTime,
           " weeks of trial time (via argument nWeeksTrialTime).\n")
    
    ## get enrollment count and times
    enr <- getEnrollment(
      rate = rateParams$enrollRate,
      enrollPeriod = enrollPeriod,
      nEnroll = nEnroll,
      nWeeks = nWeeksEnroll,
      startWeek = startWeek,
      maxEnroll= maxEnrollment )
    
    ## do a check immediately on enrollment times *if* the user
    ## specified 'nWeeksTrialTime'
    if ( !is.null( nWeeksTrialTime ) && !is.null(nEnroll) &&
           ( lastWeek <- ceiling(max(enr$times)) ) > nWeeksTrialTime ) {
      
      on.exit( cat("NOTE:  Enrollment of the", nEnroll, "participants took",
                   lastWeek, "weeks time, but only", nWeeksTrialTime, "weeks",
                   "of trial time were allowed by the user (via",
                   "argument 'nWeeksTrialTime').", 
                   "The number of enrollees is being reduced to enforce",
                   "the trial time limit specified.",
                   fill = 70, labels= "simFullEDIdata(): ") )
      
      
      ## Modify 'enr' to enforce 'nWeeksTrialTime'
      enr$times <- enr$times[ ceiling( enr$times ) <= nWeeksTrialTime ]
      enr$N <- length( enr$times )
    }
    
    ## extract enrollment count
    N <- enr$N
    
    ## generate dropout times for the 'N' participants from an exponential distribution
    dropout <- getDropout(N, rateParams$dropRate)
    
    ## generate treatment assignments
    nTrt <- length(trtAssignProb)
    
    trtAssignments <- getTreatmentAssignment(N, prob = trtAssignProb, 
                                             blockSize = blockSize )
    
    ## generate infection data for the 'N' participants
    infect <- getInfection(N, baseRate=rateParams$infecRate, 
                           relRates = infecRates, 
                           trtAssgn = trtAssignments)
    
    ## If 'nWeeksFU' or 'nWeeksTrialTime was specified, then we've got an 
    ## established bound on follow-up time.  We need to censor the dropout
    ## and infection times using that info.
    if ( !is.null( nWeeksFU ) || !is.null( nWeeksTrialTime ) )
    {
      ## amount of trial time simulated within this invocation of simulateEDI()
      if ( !is.null(nWeeksFU) ) {
        fuTime <- nWeeksFU 
      } else {
        ## else nWeeksTrialTime must have been specified
        fuTime <- (startWeek - 1 + nWeeksTrialTime) - enr$times
      }
      
      ## censor dropout and infection times by 'fuTime'
      dropout[ is.TRUE(dropout > fuTime) ] <- NA
      infect[ is.TRUE(infect > fuTime) ] <- NA
      
    } else
    {
      fuTime <- NA
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
    ## (b) the last visit occurred (during the simulated period).    ##
    ## ------------------------------------------------------------ ##
    
    ## create observed data object to fill in
    obsEDI <- data.frame( trt = trtAssignments,
                          enrollTime = enr$times,
                          dropTime = dropout,
                          infecDxTime = NA,
                          lastNegTestTime = NA,
                          futime = NA )
    
    ## get infection diagnosis dates for each non-NA infection time
    nonNA.inf <- !is.na( infect )
    
    ## Do the following only if there's at least one infection
    ## (NOTE - to reduce complexity, this section has redundant computation)
    if ( !all( is.na(infect) ) )
    {
      ## protocolVisitFunc = RSA_vstSch_3mo, which returns a vector of scheduled visit weeks assuming 3-monthly
      ## testing after the last vaccination visit
      infecDX <- getVisitWeek( infect, protocolVisits, whichVisit = "next")
      
      ## Compute the minimum of the infecDX time, dropout time and fuTime
      ## Any events occurring *strictly* after this time are censored
      minTime <- pmin( infecDX, dropout, fuTime, na.rm=TRUE)
      
      ## censor infecDX and dropout times
      infecDX[ is.TRUE(infecDX > minTime) ] <- NA
      dropout[ is.TRUE(dropout > minTime) ] <- NA
      
      ## store info in obsEDI
      obsEDI$infecDxTime <- infecDX
      obsEDI$dropTime <- dropout
    }
    
    
    ## Now that we're done making changes to 'dropout' and 'infecDxTime',
    ## we compute the amount of "follow-up time" ('futime') for each
    ## ppt.  This is equal to the 'fuTime' for ppt.s without events,
    ## and equal to the event time for participants with events.
    obsEDI$futime <- pmin( obsEDI$infecDxTime, 
                           obsEDI$dropTime, fuTime, na.rm=TRUE)
    
    
    ## Fill in the time of the last Negative HIV test.  This will be the
    ## last visit prior to 'obsTime'.  Why?  Because 'obsTime' equals
    ## infecDxTime for infecteds, dropTime for dropouts, and fuTime for
    ## everyone else, and this is what we want.
    obsEDI$lastNegTestTime <- getVisitWeek( obsEDI$futime, 
                                            protocolVisits, "previous")
    
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


## Function to determine the number (i.e. count) of the first event that
## meets all the following criteria:
##
## The arguments to the function are: 
## ----------------------------------
##   'd':
##
getFirstNonEffCnt <- 
    function(d, minCnt=75, maxCnt=Inf, lagTimes=c(26,52), lagMinCnts=c(10,2) )
{ 
  ## extract events
  events <- d[d$event==1, c("exit","entry")]

  ## order by calendar time of event (i.e. 'exit')
  events <- events[order(events$exit), ]

  ## add on a count index
  events$cnt <- 1:nrow(events)

  ## vector of infection totals at which each lag condition is met
  lagCondCnts <- mapply(FUN= function(lagtime, mincnt, dat) {
                            wLag <- ( (dat$exit - dat$entry) > lagtime)
                            dat$cnt[ wLag ][ mincnt ] },
                        lagtime=lagTimes, mincnt=lagMinCnts, 
                        MoreArgs=list(dat=events ) )

   ## The minimum of 'maxCnt' and largest of minCnt and the values that satisfy
   ## the lagMinCnts conditions
   min(maxCnt, max( minCnt, lagCondCnts ) )
}


## Function to implement "harm" monitoring on data.frame 'd' using boundaries
## specified in data.frame 'bounds'.
##
## The 'bounds' data frame must contain the columns specified by the arguments
## 'totInfecVar' and 'vaccInfecVar'.  These columns should contain a
## running total of infections in both trt groups (column 'totInfecVar'), and a 
## running total of infections in the vaccinee group ('vaccInfecVar').
## Other columns are allowed, but will be ignored.
##
## Data.frame 'd' will contain one row per infection, a 'trt' variable and a
## variable 'exit' containing the study time that the infection is observed
## (i.e. diagnosed) - the infected ppt.s "exit time" from un-infected Follow-up
do_harm_monitoring <- function(d, bounds, totInfecVar="N", vaccInfecVar="V") {
  
  ## sanity checks on input
  if ( length(d)==0 )
    stop("Argument 'd' to function 'do_harm_monitoring' has length 0.\n")
  
  if ( length(bounds)==0 )
    stop("Argument 'bounds' to function 'do_harm_monitoring' has length 0.\n")
  
  ## check that 'd' is ordered - else order it
  if ( nrow(d) > 1 & !all( diff(d$exit) >= 0 ) ) {
    d <- d[ order(d$exit), ]
  }
  
  ## get cumulative count of vaccinee infections
  vaccInfecCnt <- cumsum( d$trt > 0 ) # cumulatively counts vaccinee infections
  
  ## get count of number of infections
  totInfec <- length( vaccInfecCnt )
  
  ## restrict 'bounds' to apply to only infection totals we have
  bounds <- bounds[ bounds[, totInfecVar] <= totInfec , ]
  
  ## subset out the components of 'vaccInfecCnt' that correspond to the 
  ## infections totals at which we can stop for harm (the values in 
  ## bounds[[ totInfecVar ]]).   Compare the subsetted values to the 
  ## values in: bounds[[ vaccInfecVar ]].  If any are equal we stop.
  vaccInfecSub <- vaccInfecCnt[ bounds[[totInfecVar ]] ] # matching 'vaccInfecCnt' to 'bounds'
  harmBoundsHit <- ( vaccInfecSub >= bounds[[ vaccInfecVar ]])
  
  ## remove NA (added by Yu on 12/28/2011)
  if ( any( harmBoundsHit, na.rm=TRUE ) ) {
    
    w.hit <- which( harmBoundsHit )[1]
    
    ## 'N' is the infection count at the harm bound hit, 
    ## 'V' is the vaccine total at the harm bound hit
    N <- bounds[[totInfecVar]][ w.hit ]
    V <- bounds[[vaccInfecVar]][ w.hit ]
    
    list( isHarm = TRUE,
          stopTime = d$exit[ N ], # 'd' is ordered by 'exit' time
          stopInfectCnt = N,
          stopInfectSplit = c(Vacc = V, Plac = N - V))
    
  } else {
    list( isHarm = FALSE )
  }
}


## Function to implement "harm" monitoring on data.frame 'd' using boundaries
## specified in data.frame 'bounds'.
##
## The 'bounds' data frame must contain the columns specified by the arguments
## 'totInfecVar' and 'vaccInfecVar'.  These should columns should contain a
## running total of infections in both trt groups (column 'totInfecVar'), and a 
## running total of infections in the vaccinee group ('vaccInfecVar').
## Other columns are allowed, but will be ignored.
##
## Data.frame 'd' will contain one row per infection, a 'trt' variable and a
## variable 'exit' containing the study time that the infection is observed
## (i.e. diagnosed) - the infected ppt.s "exit time" from un-infected Follow-up
do_harm_monitoring2 <- function(d, bounds, d2, stage1 =78,  totInfecVar="N", vaccInfecVar="V") {
  
  ## sanity checks on input
  if ( length(d)==0 )
    stop("Argument 'd' to function 'do_harm_monitoring' has length 0.\n")
  
  if ( length(bounds)==0 )
    stop("Argument 'bounds' to function 'do_harm_monitoring' has length 0.\n")
  
  
  ## check that 'd' is ordered - else order it
  if ( nrow(d) > 1 & !all( diff(d$exit) >= 0 ) ) {
    d <- d[ order(d$exit), ]
  }
  
  ## get cumulative count of vaccinee infections
  vaccInfecCnt <- cumsum( d$trt > 0 )
  
  ## get count of number of infections
  totInfec <- length( vaccInfecCnt )
  
  ## restrict 'bounds' to apply to only infection totals we have
  bounds <- bounds[ bounds[, totInfecVar] <= totInfec , ]
  
  ## subset out the components of 'vaccInfecCnt' that correspond to the 
  ## infections totals at which we can stop for harm (the values in 
  ## bounds[[ totInfecVar ]]).   Compare the subsetted values to the 
  ## values in: bounds[[ vaccInfecVar ]].  If any are equal we stop.
  vaccInfecSub <- vaccInfecCnt[ bounds[[totInfecVar ]] ]
  harmBoundsHit <- ( vaccInfecSub == bounds[[ vaccInfecVar ]])
  
  ## remove NA (added by Yu on 12/28/2011)
  if ( any( harmBoundsHit, na.rm=TRUE ) ) {
    
    w.hit <- which( harmBoundsHit )[1]
    
    ## 'N' is the infection count at the harm bound hit, 
    ## 'V' is the vaccine total at the harm bound hit
    N <- bounds[[totInfecVar]][ w.hit ]
    V <- bounds[[vaccInfecVar]][ w.hit ]
    
    ## calculate the stop time at stage 1, i.e., 
    ## follow all the subjects enrolled before 'stopTime' until 18 month
    t.i = d$exit[ N ]
    D <- subset(d2, entry < t.i )  
    
    stopTimeStg1 = max(D$exit)
    
    list( isHarm = TRUE,
          stopTime = d$exit[ N ],
          stopTimeStg1 = stopTimeStg1,
          stopInfectCnt = N,
          stopInfectSplit = c(Vacc = V, Plac = N - V))
    
  } else {
    list( isHarm = FALSE )
  }
}

## Functions for harm boundary
## This function creates an object containing the values of the
## function P(n, s) which Breslow (1970, JASA) defined to be
## for 0 <= s <= n <= N), the probability of the binomial random
## walk S_n continuing to S_n = s without "absorption" in the
## rejection region (i.e. without hitting the stopping bounds).
##
## N = max number of samples
##
## S_n = sum{ x_i, i=1,..,n }, and the {x_i} are iid bernoulli(p)
##
## To calculate this, we need to specify a value of 'p', 'N', and
## a vector 'B' of length N that specifies the stopping boundares
## associated with the values of 'n' from 1:N (respectively).  So,
## B[1] will be the stopping value for n=1, B[2] the stopping value
## for n=2, ..., and B[N] the stopping value for n=N.
##
##
## The object is a list with components:
##
##   'p' - value of 'p' passed to function
##   'N' - value of 'N' passed to function
##   'Bounds' - data.frame with columns 'n' = 1:N and
##              'StoppingBound'= argument 'Bound'
##
##   'Pns' =  list of length N, each sublist containing a numeric vector.
##            Pns[[ i ]] is a numeric vector containing values of P(n,s)
##            for n=i (i.e. values of P(i, s)) for values 0<= s <= i
##
##            Note that values of s >= Bound[i] will be set to zero (see
##            defn of P(n,s) at top to see why).
##
##   'Stop' - numeric vector of length N, with Stop[i] giving the prob.
##            of hitting the stopping bound (for first time) at n=i
pNS <- function(Bound, p=.5, N=45)
{
  if( length(Bound) != N )
    stop("Length of vector 'Bound' must equal value of argument 'N'")
  
  
  ## create 'Pns' and 'Stop'
  Stop <- numeric(N)
  Pns <- vector("list", length=N)
  
  if ( is.na(Bound[1]) || Bound[1]>1 )
  {
    Pns[[ 1 ]] <- c("0" = (1-p), "1" = p)
    Stop[ 1 ] <- 0
  } else
    stop("Why are you allowing stopping when n=1!!!")
  
  for (i in 2:N)
  {
    pv <- numeric(i+1)
    names(pv) <- as.character(0:i)
    Pns[[ i ]] <- pv
    
    max.S <- min( i, Bound[i]-1, na.rm=TRUE)
    
    Pns[[i]]["0"] <- (1-p)*Pns[[i-1]]["0"]
    
    for (s in as.character(1:max.S) )
    {
      s.minus.1 <- as.character( as.numeric(s)-1)
      
      if (as.numeric(s) < i )
        Pns[[i]][ s ] <- p*Pns[[i-1]][s.minus.1] + (1-p)*Pns[[i-1]][s]
      else
        Pns[[i]][ s ] <- p*Pns[[i-1]][s.minus.1]
    }
    
    ## stopping prob. is the prob. we were one below the current bound at the
    ## last 'n', times prob. that we had another 'success'
    if ( !is.na(Bound[i]) )
      Stop[ i ] <- p*Pns[[ i-1 ]][ as.character(Bound[i]-1) ]
  }
  
  Bounds <- data.frame( n=1:N, StoppingBound=Bound )
  
  totalStopProb <- sum( Stop )
  ExpStopTime <- sum( (1:N)*(Stop) )
  
  
  list(p = p, N = N, Bounds = Bounds, Pns = Pns, Stop = Stop,
       totalStopProb = totalStopProb,  ExpStopTime = ExpStopTime)
}

####################### end of function 'pNS'#############################

### THIS USES 'nonConstBounds' framework to do constant-bounds.
### Actually, it's constant from 5 to 45, it's zero before that, so
### technically it is non-constant...


## 'x' is the total number of infections (vacc + placebo)
## 'alphaVals' is a vector of nominal (un-adjusted) p-value thresholds
##    to use in establishing cutoffs for harm-monitoring.  This vector
##    must have the same length as 'startHarmMonitor'.  The i-th value
##    of 'alphaVals' applies to the i-th interval defined by 'startHarmMonitor'
## 'startHarmMonitor' gives the endpoints of all intervals.
##    The starting point of the first interval is 1, and of the the i-th interval
##    (for i>1) is 1 + startHarmMonitor[i-1]
semiConstSpending <- function(x, alphaVals, startHarmMonitor )
{
  which.interval <- findInterval(x, c(1, startHarmMonitor),
                                 rightmost.closed=TRUE)
  alphaVals[ which.interval ]
}

getAlphaPerTest <- function(harmMonitorRange, null.p, totalAlpha=0.05) 
{
    getCumAlpha <- function( alphaPerTest, harmMonitorRange, null.p) { 
        harmBounds <- getHarmBound(
                          N = harmMonitorRange[2], 
                          per.test = alphaPerTest, 
                          harmBoundRange = harmMonitorRange,
                          null.p = null.p)
        return( harmBounds$cumStopProb[ nrow(harmBounds) ] - totalAlpha )
    }
    return( uniroot(getCumAlpha, interval = c(0,0.05), 
                    harmMonitorRange = harmMonitorRange, 
                    null.p=null.p)$root )
}

getHarmBound <- function(N,  ##Total number of infections desired for harm monitoring
                         per.test, ## value for per-test alpha level
                         harmBoundRange,
                         null.p,
                         dataDir = NULL,
                         verbose = TRUE){
  ## Note: 
  ##   'null.p' = the probability that an infection occurs in a vaccinee, under the null 
  ##              hypothesis that infection is equally likely in vaccinees and placebo
  ##              recipients.  Hence 'null.p' equals the fraction of the populations that
  ##              has received vaccine.  This would be 0.5 under a 1:1 randomization, or
  ##              11/29 under a 11:18 randomization (V:P).
  
  ## We consider the total number of infected participants to be 'N' and
  ## assume apriori and equal likelihood of infection for vaccinees as for
  ## placebos.  So the number of infected vaccinees should follow a binomial
  ## distribution with size N and probability p=null.p (where 'null.p' is 
  ## the proportion of vaccinees in the trial).  We wish to have a total
  ## probability of a "type I error" (stopping the trial for 'harm' when there
  ## is no real difference between vaccine and placebo) of .05.  Our approach
  ## will be to test for harm after each new infection, and we wish to find a
  ## fixed 'alpha' value to use for each test, so that the overall prob. of a
  ## false positive over the course of the trial is .05.
  ##
  ## To do this, we will iteratively choose a value alpha, generate a set of
  ## 'stopping bounds' corresponding to that alpha, and then estimate the
  ## overall type I error rate for those bounds.  We then go back and adjust
  ## our alpha value (up or down) depending on whether our estimated type I
  ## error is too high or too low.  Repeat until desired accuracy is obtained.
  
  bound <- NULL
  
  ## create data frame to store results in
  bounds <- data.frame(totInfec=1:N, vaccInfecBound=NA, alphaLevelBound=NA,
                       nextHigherAlphaLevel=NA, cutoff=NA )
  
  for (j in 1:nrow(bounds))
  {
    totInfec <- bounds$totInfec[j]
    
    alphaVal <- semiConstSpending( totInfec, alphaVals=c(0, per.test),
                                   startHarmMonitor = harmBoundRange)
    
    ## we don't need to do the next few steps unless alphaVal is > 0
    if (alphaVal <= 0) next
    
    ## choose the lowerBound for searching for the next cutoff value.
    ## Under our framework it will always be the same as, or higher than
    ## the cutoff from the previous (smaller) value of 'totInfec'
    if ( is.null(bound) ) {
      lowerBnd <- ceiling( null.p * totInfec )
    } else lowerBnd <- bound
    
    # startVal <- min(totInfec, bound)
    # lowerBnd <- startVal - 2
    valSeq <- totInfec:lowerBnd
    
    upperTailProbs <- cumsum( dbinom(valSeq, totInfec, null.p) )
    signif <- ( upperTailProbs <= alphaVal )
    
    ## if we have at least one significant value then do...
    if ( isTRUE(signif[1]) )
    {
      ## get "largest" (last) index for which signif == TRUE
      largest.index <- max( which( signif ) )
      
      ## define 'bound' to be the infection count corresponding to
      ## the 'largest.index' (i.e. the smallest infection count for
      ## which we have significance at per-test-level 'alphaVal'
      bound <- valSeq[ largest.index ]
      
      bounds$vaccInfecBound[ j ] <- bound
      bounds$alphaLevelBound[ j ] <- upperTailProbs[ largest.index ]
      bounds$cutoff[ j ] <- alphaVal
    }
    
  }
  out <- pNS(Bound=bounds$vaccInfecBound, p=null.p, N=N)
  
  names(out$Bounds)[ names(out$Bounds)=="StoppingBound" ] <- "Nvacc"
  boundOut <- transform(out$Bounds, Nplac= n-Nvacc, RR=round(Nvacc/(n-Nvacc),digits=2))
  boundOut <- cbind( boundOut, stopProb=round(out$Stop,4),
                     cumStopProb=round(cumsum(out$Stop),4),
                     alphaVal = bounds$cutoff )
  
  overall.alpha <- out$totalStopProb
  
  out <- pNS(Bound=bounds$vaccInfecBound, p=null.p, N=N)
  
  
  ## Add info on stopping probabilities to 'bounds' object
  bounds[, "stoppingProb"] <- out$Stop
  bounds[, "cumStoppingProb"] <- cumsum( bounds[, "stoppingProb"] )
  
  harmBounds =  boundOut
  names(harmBounds)[1:3]=c("N", "V", "P") 
  if (!is.null(dataDir)) {
      fileName <- sprintf("harmBounds_N=%d_alphaPerTest=%6.4f_pVacc=%4.2f.csv",
                          N, round(per.test, 4), round(null.p, 2) )

      write.csv(harmBounds, file.path(dataDir, fileName), row.names=FALSE)

      if (verbose) {
          cat("Potential-harm stopping boundaries saved in:\n", 
              file.path(dataDir, fileName), "\n\n") 
      }
  }
  return(harmBounds)  
}


#' Simulation of Multi-Arm Randomized Phase IIb/III  Efficacy Trials with Time-to-Event Endpoints
#'
#' \code{simTrial} generates independent time-to-event data-sets according to a user-specified trial design. The user makes assumptions about the enrollment, dropout, and infection processes in each treatment arm.
#'
#' @param N a numeric vector specifying the numbers of enrolled trial participants per treatment arm. The length of \code{N} equals the total number of treatment arms, and the first component of \code{N} represents the control arm.
#' @param aveVE a numeric vector containing, for each treatment arm in \code{N}, a time-averaged vaccine efficacy (VE), defined as the weighted average of VEs in the time intervals specified by \code{vePeriods}. If \code{VEmodel = "half"}, VE is halved in the initial interval, the full VE is applied in the second interval, and \code{aveVE} is applied thereafter. The components of \code{N} and \code{aveVE} correspond to each other.
#' @param VEmodel a character string specifying whether VE is assumed to be constant (option "\code{constant}") or have multiple levels (option "\code{half}") over time. The option "\code{half}" allows either a 2- or 3-level VE model (specified by the \code{vePeriods} vector with either 2 or 3 components). Either multi-level VE model assumes a maximal VE in the second time interval such that, when halved in the first interval, the weighted average of VE over the first two time intervals equals \code{aveVE}. Only the first character is necessary.
#' @param vePeriods a numeric vector defining start times (in weeks) of time intervals with (potentially) distinct VE levels depending on the choice of the \code{VEmodel}. If \code{VEmodel} equals "\code{half}", then \code{vePeriods} must have length 2 or 3.
#' @param enrollPeriod the final week of the enrollment period
#' @param enrollPartial the final week of the portion of the enrollment period with a reduced enrollment rate defined by \code{enrollPartialRelRate}
#' @param enrollPartialRelRate a non-negative value characterizing the fraction of the weekly enrollment rate governing enrollment from week 1 until week \code{enrollPartial}
#' @param dropoutRate a (prior) dropout probability within 1 year
#' @param infecRate a (prior) infection probability within 1 year in the control arm
#' @param fuTime a follow-up time (in weeks) of each participant
#' @param visitSchedule a numeric vector listing the visit weeks at which testing for the endpoint is conducted
#' @param missVaccProb a numeric vector with conditional probabilities of having missed a vaccination given the follow-up time exceeds \code{VEcutoffWeek} weeks. For each component, a separate per-protocol indicator is generated. Each per-protocol cohort includes subjects with (i) a non-missing vaccination, and (ii) follow-up time exceeding \code{VEcutoffWeek} weeks. If \code{NULL}, no per-protocol indicators are included.
#' @param VEcutoffWeek a time cut-off (in weeks); the follow-up time exceeding \code{VEcutoffWeek} weeks is required for inclusion in the per-protocol cohort
#' @param nTrials the number of trials to be simulated
#' @param blockSize a constant block size to be used in permuted-block randomization. The choice of \code{blockSize} requires caution to achieve the desired balance of treatment assignments within a block.
#' @param stage1 the final week of stage 1 in a two-stage trial
#' @param saveFile a character string specifying the name of the output \code{.RData} file. If \code{NULL} (default), a default file name will be used.
#' @param saveDir a character string specifying a path for the output directory. If supplied, the output is saved as an \code{.RData} file in the directory; otherwise the output is returned as a list.
#' @param verbose a logical value indicating whether information on the output directory and file name should be printed out (default is \code{TRUE})
#' @param randomSeed sets seed of the random number generator for simulation reproducibility
#'
#' @details All time variables use week as the unit of time. Month is defined as 52/12 weeks.
#'
#' The prior weekly enrollment rate is calculated based on the duration of the enrollment periods with reduced/full enrollment rates and the total number of subjects to be enrolled.
#' 
#' The weekly enrollment, dropout and infection rates used for generating trial data are sampled from specified prior distributions (the prior annual dropout and infection probabilities are specified by the user). The default choice considers non-random point-mass distributions, i.e., the prior rates directly govern the accumulation of trial data.
#'  
#' Subjects' enrollment is assumed to follow a Poisson process with a time-varying rate (the argument \code{enrollPartialRelRate} characterizes a reduced enrollment rate applied to weeks 1 through \code{enrollPartial}, i.e., full enrollment starts at week \code{enrollPartial}+1). The number of enrolled subjects is determined by the vector \code{N}.
#'  
#' Dropout times are assumed to follow an exponential distribution where the probability of a dropout within 1 week is equal to \code{dropoutRate}/52.
#'  
#' Permuted-block randomization is used for assigning treatment labels. If left unspecified by the user, an appropriate block size, no smaller than 10, will computed and used.  The function \code{getBlockSize} can be used to determine appropriate block sizes (see help(getBlockSize)).
#'  
#' Infection times are generated following the VE schedule characterized by \code{aveVE}, \code{VEmodel} and \code{vePeriods}. Independent exponential times are generated within each time period of constant VE, and their minimum specifies the right-censored infection time. Exponential rates are chosen that satisfy the user-specified requirements on the treatment- and time-period-specific probabilities of an infection within 1 week (in the control arm, the infection probability within 1 week uniformly equals \code{infecRate}/52).
#'
#'Infection diagnosis times are calculated according to the \code{visitSchedule}. The observed follow-up time is defined as the minumum of the infection diagnosis time, dropout time, and \code{fuTime}.
#'  % describe the generation of enrollment, dropout and infection rates from prior distributions to be used for data simulation (sampleRates)
#'  % describe the data generating process (each step), characterize the distributions/assumptions
#'
#' @return If \code{saveDir} is specified, the output list (named \code{trialObj}) is saved as an \code{.RData} file (the output directory path is printed); otherwise it is returned. The output object is a list with the following components:
#' \itemize{
#'   \item \code{trialData}: a list with \code{nTrials} components each of which is a \code{data.frame} with at least the variables \code{trt}, \code{entry}, \code{exit}, and \code{event} storing the treatment assignments, enrollment times, study exit times, and event indicators, respectively. The observed follow-up times can be recovered as \code{exit} - \code{entry}. Indicators of belonging to the per-protocol cohort (named \code{pp1}, \code{pp2}, etc.) are included if \code{missVaccProb} is specified.
#'   \item \code{NinfStage1}: a list whose components are numeric vectors with the numbers of \code{stage1} infections by treatment (\code{[1]} = control arm) for each simulated trial
#'   \item \code{nTrials}: the number of simulated trials
#'   \item \code{N}: the total number of enrolled trial participants
#'   \item \code{nArms}: the number of treatment arms
#'   \item \code{trtAssgnProbs}: a numeric vector containing the treatment assignment probabilities
#'   \item \code{blockSize}: the block size used for treatment assignment
#'   \item \code{fuTime}: the follow-up time (in weeks) of each participant
#'   \item \code{rates}: a list with three components: the prior weekly enrollment rate (\code{enrollment}), the prior probability of dropout within 1 week (\code{dropout}), and the prior probability of infection within 1 week (\code{infection})
#'   \item \code{enrollSchedule}: a \code{data.frame} summarizing information on enrollment periods and corresponding relative enrollment rates (relative to the weekly "base" enrollment rate). The column names are \code{start}, \code{end}, and \code{relativeRates}.
#'   \item \code{VEs}: a list with components being numeric vectors containing VE levels assumed within time periods defined by \code{vePeriods} for each active treatment arm
#'   \item \code{infecRates}: a \code{data.frame} summarizing information on time periods of distinct VE across all treatment arms. The variables \code{trt}, \code{start}, \code{end}, and \code{relRate} carry treatment assignment labels, first and last week of a time interval, and the pertaining assumed hazard ratio in the given interval.
#'   \item \code{randomSeed}: the set seed of the random number generator for simulation reproducibility
#' }
#' 
#' @examples 
#' 
#' simData <- simTrial(N=c(1000, rep(700, 2)), aveVE=seq(0, 0.4, by=0.2), 
#'                     VEmodel="half", vePeriods=c(1, 27, 79), enrollPeriod=78, 
#'                     enrollPartial=13, enrollPartialRelRate=0.5, dropoutRate=0.05, 
#'                     infecRate=0.04, fuTime=156, 
#'                     visitSchedule=c(0, (13/3)*(1:4), seq(13*6/3, 156, by=13*2/3)),
#'                     missVaccProb=c(0,0.05,0.1,0.15), VEcutoffWeek=26, nTrials=5, 
#'                     blockSize=30, stage1=78, randomSeed=300)
#'
#' ### alternatively, to save the .RData output file (no '<-' needed):
#' ###
#' ### simTrial(N=c(1400, rep(1000, 2)), aveVE=seq(0, 0.4, by=0.2), VEmodel="half", 
#' ###          vePeriods=c(1, 27, 79), enrollPeriod=78, enrollPartial=13, 
#' ###          enrollPartialRelRate=0.5, dropoutRate=0.05, infecRate=0.04, fuTime=156, 
#' ###          visitSchedule=c(0, (13/3)*(1:4), seq(13*6/3, 156, by=13*2/3)), 
#' ###          missVaccProb=c(0,0.05,0.1,0.15), VEcutoffWeek=26, nTrials=5, 
#' ###          blockSize=30, stage1=78, saveDir="./", randomSeed=300)
#'
#' @seealso \code{\link{monitorTrial}}, \code{\link{censTrial}}, and \code{\link{rankTrial}}
#' 
#' @export
simTrial <- function(N,
                    aveVE,
                    VEmodel=c("half", "constant"),
                    vePeriods,
                    enrollPeriod,
                    enrollPartial,
                    enrollPartialRelRate,
                    dropoutRate,
                    infecRate,
                    fuTime,
                    visitSchedule,
                    missVaccProb = NULL,
                    VEcutoffWeek,
                    nTrials,
                    blockSize = NULL,
                    stage1,
                    saveFile= NULL,
                    saveDir = NULL,
                    verbose = TRUE,
                    randomSeed = NULL) {

VEmodel <- match.arg(VEmodel)
  
## verify whether length of 'N' = length of 'aveVE'
if ( length(aveVE) != length(N) ) 
   stop( "Length of 'aveVE' does not match the number of treatment arms given in 'N'.\n" )

## total number of trial participants
Nppt = sum(N)

## VE for placebo arm
nullVE = aveVE[1]

## VE for vaccine arms
aveVE = aveVE[-1] 

## number of vaccine arms
nVaccArms = length(aveVE)

## total number of arms (assumes a single control arm)
nArms <- nVaccArms + 1

## create vector of treatment group names: 
trtArms <- c("C1", paste("T", 1:nVaccArms, sep=""))

## 'enrollRate' is the required weekly enrollment rate needed to enroll the expected 'Nppt' subjects 
## in 'enrollPeriod' accounting for partial enrollment rate during the initial 'enrollPartial' weeks
enrollRate <- Nppt/(enrollPartialRelRate * enrollPartial + enrollPeriod - enrollPartial)

## 'trtAssgnProbs' contains treatment assignment probabilities
trtAssgnProbs <- structure(N/Nppt, names=trtArms)

## specify block size if not provided by user
if ( is.null(blockSize) ){
  ## gets the minimum valid blockSize within the interval given by 'range' [10, Infinty)
  blockSize <- getBlockSize( N, range=c(10,Inf) )
}

## The rates in 'parSet' are "base rates" that will be 
## modified with "relative rates" for specific groups and/or time periods.

## For example, the infection rate here should be assumed infection rate
## in the population being recruited from.  Participants assigned to the
## control group and/or treatments with no effect on acquisition will
## have this base infection rate. Participants assigned to treatments
## that lower acquistion will have the base infection rate for some
## initial period of time (before the vaccination series is complete)
## and a lower rate later.

## 'parSet' contains weekly rates
parSet <- list(enrollment=enrollRate, dropout=dropoutRate/52, infection=infecRate/52)


## - 'enrollSchedule' contains information on enrollment periods and corresponding 
##    enrollment rates. 
## -'enrollSchedule' is passed as an argument to the 'simFullEDIdata' function.
## - partial enrollment within 'enrollPartial' weeks; 
## - full enrollment thereafter until the end of week 'enrollPeriod'
enrollSchedule <- data.frame(
                      start = c( 1, enrollPartial+1 ), 
                      end   = c( enrollPartial, NA ),
                      relativeRate = c( enrollPartialRelRate, 1) )

if( VEmodel!="constant" ) {
  # if the VE model is different from constant, then it is assumed to be defined by either 2 or 3 time intervals specified by 'vePeriods'
  # more than 3 time intervals are currently not supported
  if (length(vePeriods)==3){
    # this VE model assumes a VE level (variable 'vaccEff') in the second interval such that, when halved in the first interval,
    # the weighted average VE over the first two intervals equals 'aveVE';
    # VE in the third (last) interval is assumed to be 'aveVE'
    
    ## 'vaccEff' is a vector of true *full* VEs for each treatment (defined as a function of 'aveVE'), assumed in the second time interval
    # in this formula, it is also assumed that vePeriods[1]=1
    vaccEff <- aveVE * (vePeriods[3]-1)/((vePeriods[3]-1)-(vePeriods[2]-1)/2)
    
    ## If any value in 'vaccEff' is larger than 1, compute the largest possible
    ## average VE, return it as part of the error message, and stop.
    if ( any(vaccEff >= 1) ){
      
      ## which VE values are too large?
      w.too.large <-which( vaccEff >= 1 )
      
      ## compute maximum possible average VE such that 'vaccEff' is equal to 1
      ## we need to use a value that is less than 'maxEq'
      maxEq <- ((vePeriods[3]-1)-(vePeriods[2]-1)/2)/(vePeriods[3]-1)
      
      useAsMax <- floor(1000 * maxEq)/1000
      if (!(useAsMax < (maxEq - sqrt(.Machine$double.eps)))){
        useAsMax <- useAsMax - 0.001
      }
      stop("The following average vaccine efficacy(ies) is/are not attainable when",
           " using the \n", "halved VE model: ", aveVE[w.too.large], "\n", 
           "The maximum attainable VE is (approximately) ", useAsMax, "\n")
    } 
  } else if (length(vePeriods)==2){
    # this VE model assumes a VE level (variable 'vaccEff') in the second (last) interval such that, when halved in the first interval,
    # the weighted average VE over the first two intervals equals 'aveVE'
    # implication: the maximal VE must be greater in absolute value than 'aveVE' so that the average VE equals to 'aveVE'
    
    # the same formula as for the 3-level VE model except 'vePeriods[3]' is replacated with 'fuTime + 1'
    vaccEff <- aveVE * fuTime / (fuTime - (vePeriods[2] - 1) / 2)
  } else if (length(vePeriods)>3){
    stop("The argument 'vePeriods' must have length either 2 or 3 for this type of VE model.")
  }
  
} else {
  vaccEff <- aveVE
}

## - 'VEs' are true vaccine efficacies used for data generation
## - 'VEs' is a list with one component per treatment
## - each component is a vector of vaccine efficacies applied in various time
##   periods of the trial (e.g., partial-VE, full-VE, and waning-VE period)
VEs <- vector("list", nVaccArms)

for (ii in 1:nVaccArms) {
  if (VEmodel=="half"){
    if (length(vePeriods)==3){
      VEs[[ii]] = c(vaccEff[ii]/2, vaccEff[ii], aveVE[ii])      
    } else if (length(vePeriods)==2){
      VEs[[ii]] = c(vaccEff[ii]/2, vaccEff[ii])
    }
  }
  if (VEmodel=="constant"){
    VEs[[ii]] = aveVE[ii]
  }
}
names(VEs) <- paste("T", 1:length(vaccEff), sep="")

if ( length(VEs) != nVaccArms ) 
    stop( "VEs specified don't match the number of vaccine arms given in nVaccArms.\n" )

## 'infecRateTbl' contains information on relative infection rates (hazard 
##  ratios) for each treatment.  Please use "Inf" rather than NA to represent 
##  intervals that continue indefinitely.

## Each treatment must have a record starting at time 1 and must not have time
## gaps in it.  It does not need to extend to time "Inf" but typically should.
infecRateList <- vector("list", nArms)
names(infecRateList) <- c("C1", names(VEs))

infecRateList[[1]] <- data.frame( trt = "C1", start = 1, end = Inf, relRate = 1)

for (ii in 2:nArms) {
  trtName <- names(VEs)[ii-1]
  if (VEmodel=="half"){
    infecRateList[[ii]] <- data.frame( trt = trtName,  start = vePeriods,
                                       end = c(vePeriods[-1]-1, Inf), 
                                       relRate = 1 - VEs[[trtName]] )
  }
  if (VEmodel=="constant"){
    infecRateList[[ii]] <- data.frame( trt = trtName, start = vePeriods,
                                       end = Inf, relRate = 1 - VEs[[trtName]] )
  }
}
infecRateTbl <- do.call(rbind, infecRateList)

## specify prior distributions for data simulation
simParList <- parSet

## 'Constant' is a point-mass distribution (function)
simFuncList <- list(enrollment=Constant, dropout=Constant, infection=Constant) 
simPrior <- list( Function=simFuncList, params=simParList )

## set seed of random number generator
if( !is.null(randomSeed) ) {
    set.seed( randomSeed )
}

## get rates for use in generating data
## the current implementation considers constant rates, thus no need for inclusion inside the below 'for' loop
## for 'Constant' prior distributions, 'sampleRates' returns a list identical to 'parSet'
rates <- sampleRates(n=1, from=simPrior)

## create lists for storage of trial data
trialList  <- vector("list", nTrials)
infecList  <- vector("list", nTrials)
infecList2 <- vector("list", nTrials)
infecListAll <- vector("list", nTrials)
trialResult  <- vector("list", nTrials)

## 1. Generate enrollment times
## the number of enrolled subjects during a specific time interval is Poisson 
## distributed with rate = 'rate' * (end - start + 1), i.e.,
##   N <- rpois(1, lambda = rate * (end - (start-1)))
## enrollment times are uniformly distributed in the (start, end) interval, i.e.,
## runif(N, min=start-1, max=end)

## 2. Generate dropout times
## rexp(N, rate= -log( 1 - rate))

## 3. Generate infection times
## changing the 'rate' parameter to a form that gives the prob(event in interval (0,1)) = rate
## rexp(N, rate= -log( 1 - rate))

for ( i in 1:nTrials )
{
    ## generate data
    EDI.i <- simFullEDIdata(
                  rateParams = rates,
                  trtAssignProb = trtAssgnProbs,
                  blockSize = blockSize,
                  infecRates = infecRateTbl,
                  protocolVisits = visitSchedule,
                  enrollPeriod = enrollSchedule,
                  nEnroll= Nppt,
                  maxEnrollment = Nppt,
                  nWeeksFU = fuTime
              )

    ## code treatments as integers: C1=0, T1=1, T2=2...
    trtCode <- match(EDI.i$trt, trtArms) - 1


    entry <- EDI.i$enrollTime
    exit  <- entry + EDI.i$futime


    ## create flag for HIV infection
    infected <- !is.na(EDI.i$infecDxTime)

    ## make EDI.i data into a simpler form
    out <- data.frame( trt   = as.integer(trtCode),
                       entry = entry,
                       exit  = exit,
                       event = as.integer(infected)
                       )
    
    if ( !is.null(missVaccProb) ) {
      # create a set of indicators of belonging to a per-protocol cohort
      ppnames <- paste0("pp", 1:length(missVaccProb))
      for (ppIdx in 1:length(missVaccProb)) {
          ## set per-protocol indicator randomly (with prob. = 1 - missVaccProb[ppIdx])
          ## for ppt.s with follow-up beyond VEcutoffWeek
          idxNam <- ppnames[ppIdx]
          out[[idxNam]] <- 0
          FUcond <- ( EDI.i$futime > VEcutoffWeek )
          out[[idxNam]][ FUcond ] <- rbinom( sum(FUcond), 1, prob= 1 - missVaccProb[ppIdx]) 
      }
    }    
                         
    ## count number of infections by arm
    ## stg1 infections: infections occurring in first 'stage1" weeks of trial
    stg1_Infec <- infected & EDI.i$futime <= stage1

    ## add infection counts for 24-36 stage 2
    stg2_Infec <- infected & (EDI.i$futime > stage1 & EDI.i$futime <=fuTime)

    cntVec <- as.vector( table( as.factor(trtCode)[ stg1_Infec ] ) )
    cntVec2 <- as.vector( table( as.factor(trtCode)[ stg2_Infec ]))
    cntVecAll <- as.vector( table( as.factor(trtCode)[ infected ]) )
 
    infecList[[ i ]] <- cntVec
    infecList2[[ i ]] <- cntVec2
    infecListAll[[ i ]] <- cntVecAll
## store summary data into trialList
    trialList[[ i ]] <- out
}

## summary number of infections
#infectCnts <- vector("list", length=nVaccArms)
#infectCnts2 <- vector("list", length=nVaccArms)
#infectCntsAll <- vector("list", length=nVaccArms)
    
## Put everything into a "trial Object"
trialObj <- list( trialData = trialList,
                  nTrials = nTrials,
                  N = Nppt,
                  nArms = nArms,
                  trtAssgnProbs = trtAssgnProbs,
                  blockSize = blockSize,
                  fuTime = fuTime,
                  rates = parSet,
                  enrollSchedule = enrollSchedule,
                  ## time of last enrollment
                  VEs = VEs,
                  infecRates= infecRateTbl,
                  randomSeed = randomSeed,
                  NinfStage1 = infecList,
                  NinfStage2 = infecList2,
                  NinfAll = infecListAll                  
                 )

  # save trial output and information on used rates
  if (!is.null(saveDir)) {
      if (is.null(saveFile) ) 
        saveFile <- 
            paste0("simTrial_nPlac=", N[1], "_nVacc=", 
                   paste(N[-1], collapse="_"), "_aveVE=",
                   paste( round(aveVE,2), collapse="_"),
                   "_infRate=", format(infecRate, digits=3, nsmall=3),".RData" )

    save(trialObj, file=file.path(saveDir, saveFile))

    if (verbose){ 
        cat("Output saved in:\n", file.path(saveDir, saveFile), "\n\n") 
    }
  } 
  return( invisible( trialObj ) )
}

########################### End of simTrial function ###########################


## Function to determine the number (i.e. count) of the first event that
## meets all the following criteria:
##   (a) the cumulate percentage of events after "week1" having some property meets or
##       exceeds the value given in 'minPercent',  and
##   (b) the number/count of the event is at least as large as 'minCount'
##   (c) at least 2 infections after month 12 visit
##
## The arguments to the function are: 
## ----------------------------------
##   'x': an indicator vector (0/1 or FALSE/TRUE) reporting whether each
##        event has some property.  The vector should only contain info
##        on *EVENTS*. 
##
##   'minPercent': the threshold that must be met for the cumulative 
##                 percentage of "positive" entries (1s or TRUEs) in 'x'
##
##   'minCount':   the minimum event total that must be attained along
##                 with the minimum percentage threshold
##
##    week1 = 26:  first time cutoff for 'minPercent'          
##    nInfecAfterwk = 2:   number of infections required after week2
##    week2 = 52:          second time cutoff for 'nInfecAfterwk'
## 
## Example:  
##   We have defined the infection total at which we wish to perform the first
##   non-efficacy analysis for a vaccine as the smallest infection count for
##   which the percentage of infections occurring at least 6-months
##   post-enrollment reaches/exceeds 25% of infections, with at least 2 infections
##   after week 52 and the additional
##   criteria that the total must also be at least 50 (to be of sufficient 
##   size that the "large sample" normal-approximations used in the monitoring
##   criteria are reasonable).
##

getInfecCntFirstNonEff <- 
  function(  x,                 ## 'stage1' placebo and j-vaccine infections ordered by 'exit' time
             minPercent = 0,    ## min. fraction of infections that have occurred after 'week1'
             minCount = 0,      ## min. number of total infection required
             week1 = 26,        ## 
             nInfecAfterwk = 2, ## number of infections required after week2
             week2 = 52)        ## the first interim analysis can only occur after 'week2' 
  {
    ## number of infections at which first non-efficacy monitoring starts
    N1 <- NA
    
    ## time that reach 'minCnt' infections 
    nInfec <- nrow(x)
 
    futime <- x$exit - x$entry
    
    ## number of infections after 'week2'
    nPostWk2 <- sum( futime > week2 )
    
    if ( (nInfec >= minCount) && (nPostWk2 >= nInfecAfterwk) ) {
      
      ## create column containing fraction of infections post-week1
      postWk1_fract <- cumsum(futime>week1)/(1:nInfec) 
      
      ## the value of N1 is the smallest 'N' that satisfies all three conditions
      ## SIMULTANEOUSLY
      N1 <- which( (x$nInf >= minCount) & (postWk1_fract >= minPercent) &
                     ( cumsum(futime > week2) >= 2 ) )[ 1 ]
      
      ## if nothing meets all three criteria then 'N1' will have length 0
      if ( length(N1) == 0 )
        N1 <- NA
    } 
    return( N1 ) 
  }


#' Generation of Pre-Unblinded Follow-Up Data-Sets by Applying the Monitoring Outcomes
#' 
#' \code{censTrial} `correctly censors' treatment arms in data-sets generated by \code{simTrial} by including pre-unblinded follow-up data only according to the monitoring conclusions as reported by \code{monitorTrial}.
#' 
#' @param dataFile if \code{saveDir = NULL}, a list returned by \code{simTrial}; otherwise a name (character string) of an \code{.RData} file created by \code{simTrial}
#' @param monitorFile if \code{saveDir = NULL}, a list returned by \code{monitorTrial}; otherwise a name (character string) of an \code{.RData} file created by \code{monitorTrial}
#' @param stage1 the final week of stage 1 in a two-stage trial
#' @param stage2 the final week of stage 2 in a two-stage trial, i.e., the maximum follow-up time 
#' @param saveFile a character string specifying the name of the output \code{.RData} file. If \code{NULL} (default), a default file name will be used.
#' @param saveDir a character string specifying a path for both \code{dataFile} and \code{monitorFile}. If supplied, the output is also saved as an \code{.RData} file in this directory; otherwise the output is returned as a list.
#' @param verbose a logical value indicating whether information on the output directory and file name should be printed out (default is \code{TRUE})
#' 
#' @details All time variables use week as the unit of time. Month is defined as 52/12 weeks.
#' 
#' The following censoring rules are applied to each data-set generated by \code{simTrial}:
#' \itemize{
#'   \item If no vaccine arm registers efficacy or high efficacy in Stage 1, the placebo arm is censored on the date when the last vaccine arm hits the harm or non-efficacy boundary.
#'   \item If a vaccine arm hits the harm boundary, censor the arm immediately.
#'   \item If a vaccine arm hits the non-efficacy boundary, censor the arm on the earliest date of the two events: (1) the last vaccine arm hits the harm or non-efficacy boundary (if applicable); and (2) all subjects in the vaccine arm have completed the final \code{stage1} visit.
#' }
#' 
#' @return If \code{saveDir} is specified, the output list (named \code{trialListCensor}) is saved as an \code{.RData} file in \code{saveDir} (the path to \code{saveDir} is printed); otherwise it is returned. 
#' The output object is a list of length equal to the number of simulated trials, each of which is a \code{data.frame} with at least the variables \code{trt}, \code{entry}, \code{exit}, and \code{event} 
#' storing the treatment assignments, enrollment times, correctly censored study exit times, and event indicators, respectively. If available, indicators belonging to the per-protocol cohort 
#' (named \code{pp1}, \code{pp2}, etc.) are copied from the uncensored data-sets.
#' 
#' @examples
#' simData <- simTrial(N=c(1000, rep(700, 2)), aveVE=seq(0, 0.4, by=0.2), 
#'                     VEmodel="half", vePeriods=c(1, 27, 79), enrollPeriod=78, 
#'                     enrollPartial=13, enrollPartialRelRate=0.5, dropoutRate=0.05, 
#'                     infecRate=0.04, fuTime=156, 
#'                     visitSchedule=c(0, (13/3)*(1:4), seq(13*6/3, 156, by=13*2/3)),
#'                     missVaccProb=c(0,0.05,0.1,0.15), VEcutoffWeek=26, nTrials=5, 
#'                     stage1=78, randomSeed=300)
#' 
#' monitorData <- monitorTrial(dataFile=simData, stage1=78, stage2=156, 
#'                             harmMonitorRange=c(10,100), alphaPerTest=NULL, 
#'                             nonEffStartMethod="FKG", nonEffInterval=20, 
#'                             lowerVEnoneff=0, upperVEnoneff=0.4, highVE=0.7, 
#'                             stage1VE=0, lowerVEuncPower=0, alphaNoneff=0.05, 
#'                             alphaHigh=0.05, alphaStage1=0.05, 
#'                             alphaUncPower=0.05, estimand="cuminc", lagTime=26)
#'
#' censData <- censTrial(dataFile=simData, monitorFile=monitorData, stage1=78, stage2=156)
#' 
#' ### alternatively, to save the .RData output file (no '<-' needed):
#' ###
#' ### simTrial(N=c(1400, rep(1000, 2)), aveVE=seq(0, 0.4, by=0.2), VEmodel="half", 
#' ###          vePeriods=c(1, 27, 79), enrollPeriod=78, enrollPartial=13, 
#' ###          enrollPartialRelRate=0.5, dropoutRate=0.05, infecRate=0.04, fuTime=156, 
#' ###          visitSchedule=c(0, (13/3)*(1:4), seq(13*6/3, 156, by=13*2/3)), 
#' ###          missVaccProb=c(0,0.05,0.1,0.15), VEcutoffWeek=26, nTrials=30, 
#' ###          stage1=78, saveDir="./", randomSeed=300)
#' ###
#' ### monitorTrial(dataFile=
#' ###              "simTrial_nPlac=1400_nVacc=1000_1000_aveVE=0.2_0.4_infRate=0.04.RData", 
#' ###              stage1=78, stage2=156, harmMonitorRange=c(10,100), alphaPerTest=NULL, 
#' ###              nonEffStartMethod="FKG", nonEffInterval=20, lowerVEnoneff=0, 
#' ###              upperVEnoneff=0.4, highVE=0.7, stage1VE=0, lowerVEuncPower=0, 
#' ###              alphaNoneff=0.05, alphaHigh=0.05, alphaStage1=0.05, alphaUncPower=0.05, 
#' ###              estimand="cuminc", lagTime=26, saveDir="./")
#' ###
#' ### censTrial(dataFile=
#' ###          "simTrial_nPlac=1400_nVacc=1000_1000_aveVE=0.2_0.4_infRate=0.04.RData",
#' ###          monitorFile=
#' ###          "monitorTrial_nPlac=1400_nVacc=1000_1000_aveVE=0.2_0.4_infRate=0.04_cuminc.RData",
#' ###          stage1=78, stage2=156, saveDir="./")
#'  
#' @seealso \code{\link{simTrial}}, \code{\link{monitorTrial}}, and \code{\link{rankTrial}}
#' 
#' @export
censTrial <- function(dataFile,
                      monitorFile,
                      stage1,
                      stage2,
                      saveFile = NULL,
                      saveDir  = NULL,
                      verbose = TRUE){
                     
  if ( is.list(dataFile) ) {
    trialObj <- dataFile
  } else {
    if ( !is.null(saveDir) ){
      ## load in RData object (a list named 'trialObj' )
      load(file.path(saveDir, dataFile))
    } else {
      load(dataFile)
    }
  }
  
  nTrials = length(trialObj[["trialData"]])  
  nTrtArms <- as.integer( trialObj$nArm - 1 )   ## one placebo arm
  
  if ( is.list(monitorFile) ) {
    out <- monitorFile
  } else {
    if ( !is.null(saveDir) ){
      ## load in RData object (a list named 'trialObj' )
      load(file.path(saveDir, monitorFile))
    } else {
      load(monitorFile)
    }
  }
  
  ## a matrix to store bounds results from monitoring each single arm
  boundsRes = matrix(NA, ncol=nTrtArms, nrow = nTrials)
  
  ## an object to store the arm stop time when the arm stops (since trial starts), 
  stopTime = matrix(NA, ncol=nTrtArms, nrow = nTrials)
  
  ## create a list to store the trial data that all subjects have been correctly
  ## censored based on the monitoring results from all arms
  trialListCensor <- vector("list", nTrials)
  
  for (i in 1:nTrtArms) {  
    ## get the bounds 
    boundsRes [, i]= unlist( lapply( out, function(x) x[[i]]$boundHit) )
    stopTime [, i]= unlist( lapply( out, function(x) x[[i]]$stopTime) )
  }
  
  ## go through each trial
  for (i in 1:nTrials ) {
    
    ## create a list to store censored data for ith trial
    #trialCensorI = vector("list",nTrtArms+1 )
    trialCensorI = NULL
    
    ## extract data for the i-th trial
    datI <- trialObj[["trialData"]][[ i ]]
    datI$futime <- datI$exit - datI$entry
    
    ## maximum possible trial duration  
    trialStop = max(datI$exit)
    
    ## get the placebo arm
    datI.0 = subset(datI, trt == 0)  
    
    ## if none of the arms are efficacy or highefficacy (i.e. all arms are either harm or noneff),
    ## the entire trial stops when the last arm hits noneff/harm
    if (!any(boundsRes[i,]%in%c("Eff", "HighEff"))) {
      
      ## get the time when the trial stops
      trialStop = max(stopTime[i,])  
      
      ## censor the placebo arm        
      datI.0$event = datI.0$event == 1 & (datI.0$exit <= trialStop)
      datI.0$exit = pmin( datI.0$exit, trialStop)
    }
    ## if at least one arm reaches efficacy, placebo arm will continue follow up until stage 2
    ## no need of extra action, store the censored placebo
    trialCensorI = rbind(trialCensorI, datI.0)
    
    ## Now we move to censor the trt arms
    for (j in 1:nTrtArms) {
      
      ## subset jth trt arm
      datI.j <- subset(datI, trt==j )
      
      ## censor the subjects based on bounds results for all arms
      ## first, get the stop time for trial i arm j
      t.j = stopTime[i, j]
      
      ## second, check if hit harm
      if (boundsRes[i, j] %in% "Harm") {  ## if "Harm", stop
        ## censor the trt arm        
        datI.j$event = datI.j$event == 1 & (datI.j$exit <= t.j)
        datI.j$exit = pmin( datI.j$exit, t.j)
        
      } else { 
        
        if (boundsRes[i, j] %in% c("NonEffInterim", "NonEffFinal")) {  ## if hit the non efficacy bound
          ## get the stop time for this arm, 
          ## which is the end of stage 1 or when the last arm hits noneff/harm if no eff. trt arms
          endStage1 = max(datI.j$entry + stage1)
          trialStop = min (c(trialStop, endStage1))
          
          ## censor the trt arm at 'trialStop'      
          datI.j$event = datI.j$event == 1 & (datI.j$exit <= trialStop)
          datI.j$exit = pmin( datI.j$exit, trialStop)                   
          
        } else { ## hit high eff or efficacy, i will remove this later
          ## the arm will follow up to stage 2, no need of action
        }      
      }
      
      ## store the censored trt arm
      #trialCensorI [[j+1]] = datI.j
      trialCensorI = rbind(trialCensorI, datI.j)                
    }      
    trialListCensor [[i]] = trialCensorI      
  }

  ## save monitoring output
  if ( !is.null(saveDir) ) {
    if ( is.null(saveFile) ) {
        if ( is.list(monitorFile) )
          warning(
              "The output of 'censTrial' cannot be saved to a file\n",
              "You have not specified the argument 'saveFile', and a default\n",
              "filename cannot be constructed when argument 'monitorFile' is ",
              "a list.\n\n", immediate.=TRUE)
        saveFile <- paste0("trialDataCens", 
                           substr(monitorFile, 13, nchar(monitorFile)) ) 
    }
    save(trialListCensor, file = file.path(saveDir, saveFile) )

    if (verbose) { 
        cat("Trial data with correct censoring saved in:\n", 
             file.path(saveDir, saveFile), "\n\n") 
    }
  } 

  ## it should not be an "either/or" decision whether you save or output
  ## the results
  return( invisible( trialListCensor ) )
}


#' Ranking and Selection, and Head-to-Head Comparison of Individual and Pooled Treatment Arms
#' 
#' \code{rankTrial} assesses the probability of correctly selecting the winning most efficacious (individual and/or pooled) treatment arm, and assesses power to detect relative treatment efficacy in head-to-head comparisons of (individual and/or pooled) treatment arms.
#'
#' @param censFile if \code{saveDir = NULL}, a list returned by \code{censTrial}; otherwise a name (character string) of an \code{.RData} file created by \code{censTrial}
#' @param idxHighestVE an integer value identifying the treatment (vaccine) arm with the true highest VE(0--\code{stage2})
#' @param headHead a matrix (\code{ncol = 2}) of treatment arm indices for head-to-head comparisons, where the treatment with higher efficacy is listed first in each row
#' @param poolHead a matrix (\code{ncol} equals 3 or 4) of treatment arm indices for pooled-arm comparisons, where the pooled treatment with higher efficacy pooled over the first two columns is compared with the (pooled) treatment defined by columns 3 and onward. Ranking and selection of pooled arms is performed separately for each row of \code{poolHead}.
#' @param lowerVE a numeric value defining a `winning' treatment arm as one with sufficient evidence for rejecting the null hypothesis H0: VE(0--\code{stage1}) \eqn{\le} \code{lowerVE} x 100\% (typically set equal to 0)
#' @param stage1 the final week of stage 1 in a two-stage trial
#' @param stage2 the final week of stage 2 in a two-stage trial, i.e., the maximum follow-up time 
#' @param alpha one minus the nominal confidence level of the two-sided confidence interval used for testing a null hypothesis H0: VE(0--\code{stage1}) \eqn{\le} \eqn{b} x 100\% against an alternative hypothesis H1: VE(0--\code{stage1}) \eqn{>} \eqn{b} x 100\%
#' @param saveDir a character string specifying a path for \code{censFile}. If supplied, the output is also saved as an \code{.RData} file in this directory; otherwise the output is returned as a list.
#' @param verbose a logical value indicating whether information on the output directory and file name should be printed out (default is \code{TRUE})
#' 
#' @details All time variables use week as the unit of time. Month is defined as 52/12 weeks.
#' 
#' The probability of correct treatment selection is defined as the probability that the treatment arm with the highest estimated VE(0--\code{stage2}) is the one with the true highest VE(0--\code{stage2}) and, for this treatment arm, the null hypothesis H0: VE(0--\code{stage1}) \eqn{\le} \code{lowerVE} x 100\% is rejected. If \code{poolHead} is specified, the probability of correct pooled treatment selection is assessed for each set of two pooled treatment arms.
#' 
#' VE(0--\eqn{t}) is estimated as one minus the ratio of Nelson-Aalen-based cumulative incidence functions. The null hypothesis H0: VE(0--\eqn{t}) \eqn{\le} \eqn{b} x 100\% is rejected if the lower bound of the two-sided (1-\code{alpha}) x 100\% confidence interval for VE(0--\eqn{t}) lies above \eqn{b}.
#' 
#' For head-to-head individual and pooled treatment comparisons, powers to reject the null hypotheses that relative VE(0--\code{stage1}) \eqn{\le} 0\% and relative VE(0--\code{stage2}) \eqn{\le} 0\% are assessed using the aforementioned testing rule.
#' 
#' @return If \code{saveDir} is specified, the output list (named \code{out}) is saved as an \code{.RData} file in \code{saveDir} (the path to \code{saveDir} is printed); otherwise it is returned. The output object is a list with the following components:
#' \itemize{
#'   \item \code{rankSelectPw}: the probability of correct selection of the winning most efficacious individual treatment
#'   \item \code{headHeadPw}: if \code{headHead} is specified, a matrix of powers to detect relative VE(0--\code{stage1}) (column 1) and relative VE(0--\code{stage2}) (column 2) in head-to-head comparisons of individual treatment arms
#'   \item \code{poolRankSelectPw}: if \code{poolHead} is specified, a numeric vector of the probabilities of correct selection of the winning most efficacious pooled treatment for each set of pooled treatments
#'   \item \code{poolHeadPw}: if \code{poolHead} is specified, a matrix of powers to detect relative VE(0--\code{stage1}) (column 1) and relative VE(0--\code{stage2}) (column 2) in head-to-head comparisons of pooled treatment arms
#' }
#' 
#' @examples 
#' 
#' simData <- simTrial(N=c(1000, rep(700, 2)), aveVE=seq(0, 0.4, by=0.2), 
#'                     VEmodel="half", vePeriods=c(1, 27, 79), enrollPeriod=78, 
#'                     enrollPartial=13, enrollPartialRelRate=0.5, dropoutRate=0.05, 
#'                     infecRate=0.04, fuTime=156, 
#'                     visitSchedule=c(0, (13/3)*(1:4), seq(13*6/3, 156, by=13*2/3)),
#'                     missVaccProb=c(0,0.05,0.1,0.15), VEcutoffWeek=26, nTrials=5, 
#'                     stage1=78, randomSeed=300)
#'
#' monitorData <- monitorTrial(dataFile=simData, stage1=78, stage2=156, 
#'                             harmMonitorRange=c(10,100), alphaPerTest=NULL, 
#'                             nonEffStartMethod="FKG", nonEffInterval=20, 
#'                             lowerVEnoneff=0, upperVEnoneff=0.4, 
#'                             highVE=0.7, stage1VE=0, lowerVEuncPower=0, 
#'                             alphaNoneff=0.05, alphaHigh=0.05, alphaStage1=0.05, 
#'                             alphaUncPower=0.05, estimand="cuminc", lagTime=26)
#'
#' censData <- censTrial(dataFile=simData, monitorFile=monitorData, stage1=78, stage2=156)
#'                        
#' rankData <- rankTrial(censFile=censData, idxHighestVE=2, 
#'                       headHead=matrix(2:1, nrow=1, ncol=2), lowerVE=0, stage1=78, 
#'                       stage2=156, alpha=0.05)
#'
#' ### alternatively, to save the .RData output file (no '<-' needed):
#' ###
#' ### simTrial(N=c(1400, rep(1000, 2)), aveVE=seq(0, 0.4, by=0.2), VEmodel="half", 
#' ###          vePeriods=c(1, 27, 79), enrollPeriod=78, enrollPartial=13, 
#' ###          enrollPartialRelRate=0.5, dropoutRate=0.05, infecRate=0.04, fuTime=156, 
#' ###          visitSchedule=c(0, (13/3)*(1:4), seq(13*6/3, 156, by=13*2/3)), 
#' ###          missVaccProb=c(0,0.05,0.1,0.15), VEcutoffWeek=26, nTrials=30, 
#' ###          stage1=78, saveDir="./", randomSeed=300)
#' ###
#' ### monitorTrial(dataFile=
#' ###          "simTrial_nPlac=1400_nVacc=1000_1000_aveVE=0.2_0.4_infRate=0.04.RData", 
#' ###          stage1=78, stage2=156, harmMonitorRange=c(10,100), alphaPerTest=NULL, 
#' ###          nonEffStartMethod="FKG", nonEffInterval=20, 
#' ###          lowerVEnoneff=0, upperVEnoneff=0.4, highVE=0.7, stage1VE=0, 
#' ###          lowerVEuncPower=0, alphaNoneff=0.05, alphaHigh=0.05, alphaStage1=0.05, 
#' ###          alphaUncPower=0.05, estimand="cuminc", lagTime=26, saveDir="./")
#' ###
#' ### censTrial(dataFile=
#' ###  "simTrial_nPlac=1400_nVacc=1000_1000_aveVE=0.2_0.4_infRate=0.04.RData",
#' ###  monitorFile=
#' ###  "monitorTrial_nPlac=1400_nVacc=1000_1000_aveVE=0.2_0.4_infRate=0.04_cuminc.RData",
#' ###  stage1=78, stage2=156, saveDir="./")
#' ###
#' ### rankTrial(censFile=
#' ###  "trialDataCens_nPlac=1400_nVacc=1000_1000_aveVE=0.2_0.4_infRate=0.04_cuminc.RData",
#' ###  idxHighestVE=2, headHead=matrix(2:1, nrow=1, ncol=2), lowerVE=0, stage1=78, 
#' ###  stage2=156, alpha=0.05, saveDir="./")
#' 
#' @seealso \code{\link{simTrial}}, \code{\link{monitorTrial}}, and \code{\link{censTrial}}
#' 
#' @export
rankTrial <- function(censFile,
                      idxHighestVE,
                      headHead=NULL,
                      poolHead=NULL,
                      lowerVE,
                      stage1,
                      stage2,
                      alpha,
                      saveDir = NULL,
                      verbose = TRUE){                 
  
  if (!is.null(saveDir)){
    ## load the trial data in RData object (a list named 'trialObj' )
    load(file.path(saveDir, censFile))
  } else {
    trialListCensor <- censFile
    rm(censFile)
  }  
  
  nTrials = length(trialListCensor)  
  nTrtArms <- length(unique(trialListCensor[[1]]$trt)) - 1    ## one placebo arm
  
  if (!is.null(headHead)){
    if (NCOL(headHead)!=2){ stop("Number of columns in headHead must equal to 2.'\n") }
    
    ## a matrix to store power for head to head comparison of single arm
    ## column 1 = power for VE(0-18), column 2 = power for VE(0-36)
    headHeadPw = matrix(0, nrow=NROW(headHead), ncol=2)    
  }  
  
  if (!is.null(poolHead)){
    if (!(NCOL(poolHead) %in% 3:4)){ stop("Number of columns in poolHead must equal to 3 or 4.'\n") }
    
    ## a matrix to store power for comparison of pooled arms
    ## column 1 = power for VE(0-18), column 2 = power for VE(0-36)
    poolHeadPw = matrix(0, nrow=NROW(poolHead), ncol=2)    
    
    # power to correctly identify the best pooled vaccine arm
    poolRankSelectPw <- numeric(NROW(poolHead))
  }  
  
  # power to correctly identify the best vaccine arm
  rankSelectPw <- 0
  
  ## go through each trial
  for (i in 1:nTrials ) {
    ## extract data for the i-th trial
    datI <- trialListCensor[[i]] 
    
    ## a matrix to store estimated VE for each trt arm 
    ## column 1 = power for VE(0-18), column 2 = power for VE(0-36)
    VE.I = matrix(0, nrow=nTrtArms, ncol=2)
    
    ## cum hazard-based Wald test results for each trt arm
    ## column 1 = power for VE(0-18), column 2 = power for VE(0-36)
    test.I = matrix(0, nrow=nTrtArms, ncol=2)
    
    ## now calculate power for rank and select
    for (j in 1:nTrtArms){
      ## extract data for the j-th trt arm & placebo
      testData = subset(datI, trt %in% c(0, j))
      
      ## power for stage 1
      bnd = testVE(datI=testData, lowerVE=lowerVE, stage =stage1, alpha=alpha)
      VE.I[j, 1] = bnd$VE
      test.I[j, 1] = bnd$bound
      
      ## power for stage 2
      bnd = testVE(datI=testData, lowerVE=lowerVE, stage =stage2, alpha=alpha)
      VE.I[j, 2] = bnd$VE
      test.I[j, 2] = bnd$bound
    }
    
    # identify the vaccine arm with the highest estimated VE(0-36)
    bestInd = which(VE.I[,2] == max(VE.I[,2]))
      
    # to estimate the probability that the vaccine arm with the highest estimated VE(0-36) is the one with the
    # true highest VE(0-36), AND for that regimen the hypothesis H0: VE(0-18)<=0% was rejected    
    if (bestInd==idxHighestVE && test.I[bestInd, 1]=="Eff"){ rankSelectPw = rankSelectPw + 1 }    
    
    if (!is.null(headHead)){  
      # head-to-head comparison of individual vaccine arms at 18 and 36 months
      pw = headTestVE(datI, headHead, stage1, stage2, alpha)
      headHeadPw = headHeadPw + pw  # counts of detected relative efficacy for each comparison            
    }    
    
    if (!is.null(poolHead)) {
      # head-to-head comparison of pooled vaccine arms at 18 and 36 months
      pw = headTestVE(datI, poolHead, stage1, stage2, alpha)
      poolHeadPw = poolHeadPw + pw   
      
      # ranking and selection for pooled vaccine arms
      for (hi in 1:NROW(poolHead)) {
        arm1 = poolHead[hi,1:2]
        arm2 = poolHead[hi,3:NCOL(poolHead)]
        
        ## change trt index to combine arms listed in "arm1" and "arm2"
        ## and store the new data in datI.p
        datI$trt[datI$trt %in% arm1] = arm1[1]
        
        ## if arm2 is placebo, then no changes
        datI$trt[datI$trt %in% arm2] = arm2[1]
        
        ## now update 'bestVE' to arm1[1]
        bestVE = arm1[1]
        
        ## a matrix to store estimated VE for each trt arm 
        ## column 1 = power for VE(0-18), column 2 = power for VE(0-36)
        
        ## make it very small and impossible to be selected
        VE.I = matrix(-10, nrow=nTrtArms, ncol=2)
        test.I = matrix("empty", nrow=nTrtArms, ncol=2)
        
        ## get index for combined arms
        indTrt = sort(unique(datI$trt))
        indTrt = indTrt[-1]    ## remove placebo
        
        for (j in indTrt) {
          ## extract data for the j-th trt arm & placebo
          testData = subset(datI, trt %in% c(0, j)) 
          
          ## power for stage 1
          bnd = testVE(datI=testData, lowerVE=lowerVE, stage=stage1, alpha=alpha)
          VE.I[j, 1] = bnd$VE
          test.I[j, 1] = bnd$bound
          
          ## power for stage 2
          bnd = testVE(datI=testData, lowerVE=lowerVE, stage=stage2, alpha=alpha)
          VE.I[j, 2] = bnd$VE
          test.I[j, 2] = bnd$bound
        }
        
        # identify the pooled vaccine arm with the highest estimated VE(0-36)
        bestInd = which(VE.I[,2]==max(VE.I[,2]))
        
        # to estimate the probability that the pooled vaccine arm with the highest estimated VE(0-36) is the one 
        # with the true highest VE(0-36), AND for that regimen the hypothesis H0: VE(0-18)<=0% was rejected    
        if (bestInd==bestVE && test.I[bestInd, 1]=="Eff"){ poolRankSelectPw[hi] = poolRankSelectPw[hi] + 1 }
      }
    }
  }
  
  out <- list(rankSelectPw = rankSelectPw/nTrials)
  if (!is.null(headHead)){ out$headHeadPw <- headHeadPw/nTrials } 
  if (!is.null(poolHead)){ 
    out$poolRankSelectPw <- poolRankSelectPw/nTrials
    out$poolHeadPw <- poolHeadPw/nTrials    
  }
  
  if (!is.null(saveDir)){
    saveFile <- paste0("rankSelectPwr",substr(censFile, 14, nchar(censFile)))
    save(out, file = file.path(saveDir, saveFile) )
    if (verbose){ cat("Output saved in:\n", file.path(saveDir, saveFile), "\n\n") }
  } else {
    return(out)
  }  
}

testVE <- function(datI, lowerVE, stage, alpha){
  upperFR <- 1-lowerVE
  ## variables for cumulative hazard-based Wald test, censor at 'stage' 
  ## which can be stage 1 or 2 
  datI$futime <- datI$exit - datI$entry
  datI$event <- datI$event == 1 & ( datI$futime <= stage)
  datI$exit  <- pmin( datI$exit, datI$entry + stage )
  datI$futime <- datI$exit - datI$entry
  
  ## convert 'trt' to indicator variable before calculating VE 
  ## (i.e. convert the non-zero values to 1)
  datI$trt <- as.integer(datI$trt > 0 )
  
  # if a vaccine regimen is highly efficacious, there will be no events observed in this arm and the Nelson-Aalen estimator will be incalculable
  if (sum(datI[datI$trt==1,"event"])==0){
    FR <- 0
    bound <- "Eff"
  } else {
    KM <- survfit(Surv(futime, event) ~ trt, data=datI, error="greenwood")
    KM.sum <- summary(KM)
    # Nelson-Aalen estimates
    na.0 <- cumsum(KM.sum$n.event[KM.sum$strata=="trt=0"]/KM.sum$n.risk[KM.sum$strata=="trt=0"])
    varna.0 <- cumsum(KM.sum$n.event[KM.sum$strata=="trt=0"]/KM.sum$n.risk[KM.sum$strata=="trt=0"]^2)
    na.0 <- na.0[length(na.0)]
    varna.0 <- varna.0[length(varna.0)]        
    na.1 <- cumsum(KM.sum$n.event[KM.sum$strata=="trt=1"]/KM.sum$n.risk[KM.sum$strata=="trt=1"])
    varna.1 <- cumsum(KM.sum$n.event[KM.sum$strata=="trt=1"]/KM.sum$n.risk[KM.sum$strata=="trt=1"]^2)
    na.1 <- na.1[length(na.1)]
    varna.1 <- varna.1[length(varna.1)]
    
    # survival estimates
    S.0 <- exp(-na.0)
    varS.0 <- ifelse(is.na(varna.0), NA, exp(-2*na.0) * varna.0)
    S.1 <- exp(-na.1)
    varS.1 <- ifelse(is.na(varna.1), NA, exp(-2*na.1) * varna.1)
    
    # cumulative incidence ratio
    F.0 <- 1 - S.0
    F.1 <- 1 - S.1
    FR <- F.1/F.0
    varlogFR <- ifelse(is.na(varS.0) | is.na(varS.1), NA, varS.1/(F.1^2) + varS.0/(F.0^2))
    
    FRci.up <- ifelse(is.na(varlogFR), NA, exp(log(FR)+qnorm(1-alpha/2)*sqrt(varlogFR)))
    bound <- ifelse(FRci.up < upperFR, "Eff", "nonEff")  
  }
  
  return(list(bound=bound, VE=1-FR))
}

headTestVE <- function(datI, headHeadInd, stage1, stage2, alpha){
  ## a matrix to store power for head to head comparison 
  ## column 1 = power for VE(0-18), column 2 = power for VE(0-36)
  headHeadPw = matrix(0, nrow=NROW(headHeadInd), ncol=2)
  
  for (h in 1:NROW(headHeadInd)) {
    if(NCOL(headHeadInd)==2) {  ## 2 columns
      arm1 = headHeadInd[h,1]
      arm2 = headHeadInd[h,2]
    } else {                    ## 3 or 4 column
      arm1 = headHeadInd[h,1:2]
      arm2 = headHeadInd[h,3:NCOL(headHeadInd)]
    }
    
    testDataRel <- subset(datI, trt %in% c(arm1, arm2))
    ## need to change the index in arm2 to 0
    testDataRel$trt [testDataRel$trt %in% arm2] = 0    
    
    testData <- subset(datI, trt %in% c(0, arm1))
        
    # head-to-head power = the 1-sided head-to-head test rejects; AND the superior arm passes 
    # Stage 1 with VE(0-stage1) > 0%
    ## relative VE(0-stage1)    
    bndRel <- testVE(datI=testDataRel, lowerVE=0, stage=stage1, alpha=alpha)
    ## VE(0-stage1) of the superior treatment
    bnd <- testVE(datI=testData, lowerVE=0, stage=stage1, alpha=alpha)    
    if (bndRel$bound =="Eff" & bnd$bound=="Eff"){ headHeadPw[h, 1] = 1 }
    
    ## power for stage 2
    bndRel <- testVE(datI=testDataRel, lowerVE=0, stage=stage2, alpha=alpha)
    bnd <- testVE(datI=testData, lowerVE=0, stage=stage2, alpha=alpha)    
    if (bndRel$bound =="Eff" & bnd$bound=="Eff"){ headHeadPw[h, 2] = 1 }
  }
  return(headHeadPw)
}


#########
buildBounds = function(nInfec, highEffBounds) {
   
   ## high efficacy : approxiamtion 
   idx = which(highEffBounds$N==nInfec)
   
   ## get number of row
   nBounds = nrow (highEffBounds)
   if (length(idx)==0) {  ## no bound for nInfec, add it
      highEffBounds  = rbind(highEffBounds, c(nInfec, highEffBounds[nBounds,2]))   ## approximation  
      
      ## order the bound by infections
      highEffBounds = highEffBounds[order(highEffBounds$N), ]
      
      idx = which(highEffBounds$N==nInfec)
      if(idx>1 && (idx+1)<=nrow (highEffBounds))       ## approximate the bound by averaging over 2 adjacent twos
         highEffBounds[idx,2] = mean(c(highEffBounds[idx-1,2], highEffBounds[idx+1,2]))
    }
    highEffBounds = subset(highEffBounds, N<=nInfec)
 
   list(highEffBounds=highEffBounds)            
}


#' Unconditional Power to Detect Positive Treatment Efficacy in a Per-Protocol Cohort
#' 
#' \code{VEpowerPP} computes unconditional power to detect positive treatment (vaccine) efficacy in per-protocol cohorts identified in \code{simTrial}-generated data-sets.
#'                  
#' @param dataList if \code{saveDir = NULL}, a list of objects (lists) returned by \code{censTrial}; otherwise a list of \code{.RData} file names (character strings) generated by \code{censTrial}
#' @param lowerVEuncPower a numeric value specifying a one-sided null hypothesis H0: VE(\code{VEcutoffWeek}--\code{stage1}) \eqn{\le} \code{lowerVEuncPower} x 100\%. Unconditional power (i.e., accounting for sequential monitoring) to reject H0 in the per-protocol cohort is calculated, where the rejection region is defined by the lower bound of the two-sided (1-\code{alphaUncPower}) x 100\% confidence interval for VE(\code{VEcutoffWeek}--\code{stage1}) being above \code{lowerVEuncPower} (typically a number in the 0--0.5 range).
#' @param alphaUncPower one minus the nominal confidence level of the two-sided confidence interval used to test the one-sided null hypothesis H0: VE(\code{VEcutoffWeek}--\code{stage1}) \eqn{\le} \code{lowerVEuncPower} x 100\% against the alternative hypothesis H1: VE(\code{VEcutoffWeek}--\code{stage1}) \eqn{>} \code{lowerVEuncPower} x 100\%.
#' @param VEcutoffWeek a cut-off time (in weeks). Only subjects with the follow-up time exceeding \code{VEcutoffWeek} are included in the per-protocol cohort.
#' @param stage1 the final week of stage 1 in a two-stage trial
#' @param outName a character string specifying the output \code{.RData} file name. If \code{outName = NULL} but \code{saveDir} is specified, the output file is named \code{VEpwPP.RData}.
#' @param saveDir a character string specifying a path for the output directory. If supplied, the output is saved as an \code{.RData} file named \code{outName} in the directory; otherwise the output is returned as a list.
#' @param verbose a logical value indicating whether information on the output directory and file name should be printed out (default is \code{TRUE})
#' 
#' @details All time variables use week as the unit of time. Month is defined as 52/12 weeks.
#' 
#' A per-protocol cohort indicator is assumed to be included in the \code{simTrial}-generated data-sets, which is ensured by specifying the \code{missVaccProb} argument in \code{simTrial}.
#' 
#' VE(\code{VEcutoffWeek}--\code{stage1}) is estimated as one minus the ratio of Nelson-Aalen-based cumulative incidence functions. \code{VEpowerPP} computes power to reject the null hypothesis H0: VE(\code{VEcutoffWeek}--\code{stage1}) \eqn{\le} \code{lowerVEuncPower} x 100\%. H0 is rejected if the lower bound of the two-sided (1-\code{alphaUncPower}) x 100\% confidence interval for VE(\code{VEcutoffWeek}--\code{stage1}) lies above \code{lowerVEuncPower}.
#'
#' @return If \code{saveDir} is specified, the output list (named \code{pwList}) is saved as an \code{.RData} file named \code{outName} (or \code{VEpwPP.RData} if left unspecified); otherwise the output list is returned. The output object is a list (of equal length as \code{dataList}) of lists with the following components:
#' \itemize{
#'   \item \code{VE}: a numeric vector of VE(\code{VEcutoffWeek}--\code{stage1}) estimates for each missing vaccination probability in \code{missVaccProb} of \code{simTrial}
#'   \item \code{VEpwPP}: a numeric vector of powers to reject the null hypothesis H0: VE(\code{VEcutoffWeek}--\code{stage1}) \eqn{\le} \code{lowerVEuncPower} x 100\% for each missing vaccination probability in \code{missVaccProb} of \code{simTrial}
#' }
#' 
#' @examples 
#' simData <- simTrial(N=rep(1000, 2), aveVE=c(0, 0.4), VEmodel="half", 
#'                     vePeriods=c(1, 27, 79), enrollPeriod=78, 
#'                     enrollPartial=13, enrollPartialRelRate=0.5, dropoutRate=0.05, 
#'                     infecRate=0.04, fuTime=156, 
#'                     visitSchedule=c(0, (13/3)*(1:4), seq(13*6/3, 156, by=13*2/3)),
#'                     missVaccProb=c(0,0.05,0.1,0.15), VEcutoffWeek=26, nTrials=5, 
#'                     stage1=78, randomSeed=300)
#'
#' monitorData <- monitorTrial(dataFile=simData, stage1=78, stage2=156, 
#'                             harmMonitorRange=c(10,100), alphaPerTest=NULL, 
#'                             nonEffStartMethod="FKG", nonEffInterval=20, 
#'                             lowerVEnoneff=0, upperVEnoneff=0.4, 
#'                             highVE=0.7, stage1VE=0, lowerVEuncPower=0, 
#'                             alphaNoneff=0.05, alphaHigh=0.05, alphaStage1=0.05, 
#'                             alphaUncPower=0.05, estimand="cuminc", lagTime=26)
#'
#' censData <- censTrial(dataFile=simData, monitorFile=monitorData, stage1=78, stage2=156)
#'
#' VEpwPP <- VEpowerPP(dataList=list(censData), lowerVEuncPower=0, alphaUncPower=0.05,
#'                     VEcutoffWeek=26, stage1=78)
#'
#' ### alternatively, to save the .RData output file (no '<-' needed):
#' ###
#' ### simTrial(N=rep(1000, 2), aveVE=c(0, 0.4), VEmodel="half", 
#' ###          vePeriods=c(1, 27, 79), enrollPeriod=78, enrollPartial=13, 
#' ###          enrollPartialRelRate=0.5, dropoutRate=0.05, infecRate=0.04, fuTime=156, 
#' ###          visitSchedule=c(0, (13/3)*(1:4), seq(13*6/3, 156, by=13*2/3)), 
#' ###          missVaccProb=c(0,0.05,0.1,0.15), VEcutoffWeek=26, nTrials=5, 
#' ###          stage1=78, saveDir="./", randomSeed=300)
#' ###
#' ### monitorTrial(dataFile=
#' ###          "simTrial_nPlac=1000_nVacc=1000_aveVE=0.4_infRate=0.04.RData", 
#' ###          stage1=78, stage2=156, harmMonitorRange=c(10,100), alphaPerTest=NULL, 
#' ###          nonEffStartMethod="FKG", nonEffInterval=20, 
#' ###          lowerVEnoneff=0, upperVEnoneff=0.4, highVE=0.7, stage1VE=0, 
#' ###          lowerVEuncPower=0, alphaNoneff=0.05, alphaHigh=0.05, alphaStage1=0.05, 
#' ###          alphaUncPower=0.05, estimand="cuminc", lagTime=26, saveDir="./")
#' ###
#' ### censTrial(dataFile=
#' ###  "simTrial_nPlac=1000_nVacc=1000_aveVE=0.4_infRate=0.04.RData",
#' ###  monitorFile=
#' ###  "monitorTrial_nPlac=1000_nVacc=1000_aveVE=0.4_infRate=0.04_cuminc.RData",
#' ###  stage1=78, stage2=156, saveDir="./")
#' ###
#' ### VEpowerPP(dataList=
#' ###  list("trialDataCens_nPlac=1000_nVacc=1000_aveVE=0.4_infRate=0.04_cuminc.RData"),
#' ###  lowerVEuncPower=0, alphaUncPower=0.05, VEcutoffWeek=26, stage1=78, saveDir="./")
#' 
#' @seealso \code{\link{simTrial}}
#' 
#' @export
VEpowerPP <- function( dataList, 
                       lowerVEuncPower, 
                       alphaUncPower, 
                       VEcutoffWeek,
                       stage1, 
                       outName = NULL, 
                       saveDir = NULL, 
                       verbose = TRUE){

  upperFRuncPower <- 1-lowerVEuncPower
  # output list (for each censTrial output object) of lists with components 'VE' and 'VEpwPP'
  pwList <- as.list(NULL)
  for (k in 1:length(dataList)){
    if (!is.null(saveDir)){
      # assumes 'dataList[[k]]' is a character string
      load(file.path(saveDir, dataList[[k]]))
      dataList[[k]] <- trialListCensor
      rm(trialListCensor)
    }
    
    nTrials <- length(dataList[[k]])
    nPPcohorts <- length(grep("pp", colnames(dataList[[k]][[1]])))
    if (nPPcohorts==0){ stop("Missing per-protocol cohort indicator in the data-set.") }
    ppnames <- paste0("pp", 1:nPPcohorts)
    
    VE <- matrix(NA, nrow=nTrials, ncol=nPPcohorts)
    VEpwPP <- matrix(NA, nrow=nTrials, ncol=nPPcohorts)
    for (i in 1:nTrials){
      dataI <- dataList[[k]][[i]]
      dataI$futime <- dataI$exit - dataI$entry
      # count only infections between VEcutoffWeek and the end of Stage 1
      dataI$event <- dataI$event==1 & dataI$futime>VEcutoffWeek & dataI$futime<=stage1
      # censor everyone at the end of Stage 1
      dataI$futime <- pmin(dataI$futime, stage1)

      # AN ALTERNATIVE WAY TO COMPUTE PER-PROTOCOL VE      
#       # restrict to subjects with follow-up time exceeding 'VEcutoffWeek' weeks (per-protocol criterion 1)
#       dataI <- subset(dataI, futime > VEcutoffWeek)
#       # censor all subjects at the Week 'stage1' visit
#       dataI$event <- dataI$event == 1 & (dataI$futime <= stage1)
#       dataI$futime <- pmin(dataI$futime, stage1)
#       # shift the time origin to the Week 'VEcutoffWeek' visit
#       dataI$futime <- dataI$futime - VEcutoffWeek      
      
      for (j in 1:nPPcohorts){
        # restrict to subjects with non-missing vaccinations (per-protocol criterion 2)
        dataIj <- dataI[dataI[[ppnames[j]]]==1,]
        # now we have the final per-protocol data-set at the end of 'stage 1'
        # get the group sizes in the PP cohort
        nAll <- table(dataIj$trt)
        if (length(nAll)<2){
          next
        } else {
          # get the infection counts in the PP cohort
          nInf <- table(dataIj$trt[dataIj$event==1])
          if (length(nInf)<2){
            next
          } else {
            KM <- survfit(Surv(futime, event) ~ trt, data=dataIj, error="greenwood")
            KM.sum <- summary(KM)
            
            # Nelson-Aalen estimates
            na.0 <- cumsum(KM.sum$n.event[KM.sum$strata=="trt=0"]/KM.sum$n.risk[KM.sum$strata=="trt=0"])
            varna.0 <- cumsum(KM.sum$n.event[KM.sum$strata=="trt=0"]/KM.sum$n.risk[KM.sum$strata=="trt=0"]^2)
            na.0 <- na.0[length(na.0)]
            varna.0 <- varna.0[length(varna.0)]
            
            na.1 <- cumsum(KM.sum$n.event[KM.sum$strata=="trt=1"]/KM.sum$n.risk[KM.sum$strata=="trt=1"])
            varna.1 <- cumsum(KM.sum$n.event[KM.sum$strata=="trt=1"]/KM.sum$n.risk[KM.sum$strata=="trt=1"]^2)
            na.1 <- na.1[length(na.1)]
            varna.1 <- varna.1[length(varna.1)]
            
            # survival estimates
            S.0 <- exp(-na.0)
            varS.0 <- ifelse(is.na(varna.0), NA, exp(-2*na.0) * varna.0)
            S.1 <- exp(-na.1)
            varS.1 <- ifelse(is.na(varna.1), NA, exp(-2*na.1) * varna.1)
            
            # cumulative incidence ratio
            F.0 <- 1 - S.0
            F.1 <- 1 - S.1
            FR.i <- F.1/F.0
            varlogFR <- ifelse(is.na(varS.0) | is.na(varS.1), NA, varS.1/(F.1^2) + varS.0/(F.0^2))
            
            # VE(6.5-stage1) estimate
            VE[i,j] <- 1 - FR.i
            
            # 1-sided test of null hypothesis VE(6.5-stage1)<=lowerVEuncPower
            FRci.up <- ifelse(is.na(varlogFR), NA, exp(log(FR.i)+qnorm(1-alphaUncPower/2)*sqrt(varlogFR)))
            # if the upper bound of the CI for the cumulative incidence ratio lies below 'upperFRuncPower', reject H0
            VEpwPP[i,j] <- ifelse(FRci.up < upperFRuncPower, 1, 0)            
          }          
        }        
      }
    }
    VE <- apply(VE, 2, mean, na.rm=TRUE)
    VEpwPP <- apply(VEpwPP, 2, mean, na.rm=TRUE)
    pwList[[k]] <- list(VE=VE, VEpwPP=VEpwPP)    
  }
  if (!is.null(saveDir)){
    saveFile <- ifelse(is.null(outName), "VEpwPP.RData", outName)
    save(pwList, file=file.path(saveDir, saveFile))
    if (verbose){ cat("Output saved in:\n", file.path(saveDir, saveFile), "\n\n") }
  } else {
    return(pwList)
  }  
}


## compute CIR on given dataset, at given follow-up times
## If you want the CIR at particular 'calendar' times, you should first use
## censorTrials to create the appropriate datasets and then apply cumIncRatio
## to those censored trials.
##
cumIncRatio <- function( d, times="last", arms=c(1,0), alphaLevel=0.05, 
                         randFraction = 0.5, ..., SIMPLIFY=TRUE) 
{
    ##     d: a data.frame or list of data.frames 
    ##
    ## times: vector that can take values:
    ##          "last" (the default)- This directs that the CIR be provided for 
    ##               the latest event time only.
    ##          "all"  - directs that the CIR be provided for all event times
    ##           (NOT IMPLEMENTED ) numeric vector - interpreted as the follow-up times for which
    ##               the CIR should be provided
    ##
    ##  arms   numeric vector of length 2.   Indicates which treatment arms to 
    ##         use in computing the CIR.  The values should be the number (or
    ##         names, if assigned), for the treatments.  The first value of arms
    ##         specifies the treatment to be used in the numerator of the ratio,
    ##         and the 2nd value the treatment to be compared to (i.e. in the 
    ##         treatment used in the denominator).
    ##
    ## alphaLevel: the two-sided alpha-level used to construct confidence 
    ##             intervals for the estimate
    ##
    ## randFraction: the fraction of randomizations going to the particular
    ##               treatment group, when only that treatment group and the
    ##               placebo group are considered
    ##
    ##  ...    Not used at present.  Including this argument force exact 
    ##         specification of argument 'SIMPLIFY' if user desires to change
    ##         it (partial matching doesn't work after "..." arguments)
    ##       
    ##  SIMPLIFY  Logical - Should the return structure be simplified from a 
    ##            list of data.frames to a single data.frame for cases where
    ##            either:
    ##              - a single data set was specified by argument 'd',  or
    ##              - a single timepoint was specified for each dataset 
    ##
    ##  NOTE!! If there are no infections in one of the two groups at any of the
    ##         time-points, an adhoc procedure will be used to compute the CI for
    ##         the VE estimate, which is based on writing the VE as a function of
    ##         the ratio of infections in the two groups, treating the infection
    ##         counts as binomial, and utilizing an exact CI for the binomial 'p'
    ##         (see notes for function 'VEci_binom')
   
    ## force evaluation of objects so that they can be successfully passed
    ## through to other functions
    force( d )
    force( alphaLevel )
    force( randFraction )

    ## survival model formula to use
    survForm <- Surv(exit-entry, event) ~ trt

    ## is 'd' is a single data.frame, convert it into a list 
    if (is.data.frame(d))  
      d <- list(d)


    cir <- function(d.i, times, formula, arms, alpha, 
                    randFraction, simplify=TRUE) {

        force(d.i)
        #force(alpha)

        SurvSumm <- summary( survfit(formula, data=d.i) ) 
        SurvSumm <- as.data.frame( SurvSumm[c("strata","time","n.event","n.risk")] )
             
        ## if there are no events in either group, then punt
        ## fill in return object with NAs and return
        if ( nrow(SurvSumm)==0 ) {
            return( data.frame(
                        evalTime  = NA,
                        nEvents   =  0,
                        nEvents.1 =  0,
                        nEvents.2 =  0,
                        FR        = NA,
                        varlogFR  = NA,
                        FR_loCI   = NA,
                        FR_upCI   = NA,
                        VE        = NA,
                        VE_loCI   = NA,
                        VE_upCI   = NA ) )
        }

        strataA <- ( SurvSumm$strata == paste0("trt=", arms[1]) )
        timesA  <- SurvSumm$time[strataA]
        trtA    <- SurvSumm[ strataA , c("n.event", "n.risk")]

        strataB <- ( SurvSumm$strata == paste0("trt=", arms[2]) )
        timesB  <- SurvSumm$time[strataB]
        trtB    <- SurvSumm[ strataB , c("n.event", "n.risk")]

        NelsAal.A <- cumsum( trtA$n.event / trtA$n.risk )
        NelsAal.B <- cumsum( trtB$n.event / trtB$n.risk )

        varNelAal.A <- cumsum( trtA$n.event / trtA$n.risk^2 )
        varNelAal.B <- cumsum( trtB$n.event / trtB$n.risk^2 )

        F.A <- 1 - exp(-NelsAal.A)
        F.B <- 1 - exp(-NelsAal.B)

        varF.A <- exp(-2*NelsAal.A)*varNelAal.A
        varF.B <- exp(-2*NelsAal.B)*varNelAal.B

        ## Determine times at which to compute CIR, if not provided
        if ( times %in% c("last", "all")  ) {

            EventTimes <- sort( unique( SurvSumm$time ) )
            nTimes <- length(EventTimes)

            if (times == "all") {
                whichTimes <- 1:nTimes
            } else {
                whichTimes <- nTimes
            }
            times <- EventTimes[ whichTimes ]
        } 

        ## The 'cut' function (with labels=FALSE) returns the number of the 
        ## interval (defined by 'breaks') into which each of the 'times' falls
        ##
        ## For the breaks, we need to add "0" as a lower bound in case someone
        ## specifies a time before the first event.  This extra interval throws
        ## off the previously 1:1 correspondence with the F.A, varNelAal.A, etc.
        ## Therefore, before using the objects created using 'cut' for 
        ## subsetting we'll prepend the vectors being subsetted (see usage in
        ## code blocks following the one immediately below)
        interval.A <- cut( times, breaks=c(0, timesA, Inf), 
                           include.lowest=TRUE, right=FALSE, labels=FALSE) 
        interval.B <- cut( times, breaks=c(0, timesB, Inf), 
                           include.lowest=TRUE, right=FALSE, labels=FALSE) 

        ## values of F.A and F.B and the times provided by 'times'
        F.A_times <- c(0, F.A)[interval.A]
        F.B_times <- c(0, F.B)[interval.B] 

        FR_times    <- F.A_times/F.B_times
        logFR_times <- log( FR_times )

        ## pull out values needed for computing variance of logFR_times
        varF.A_times <- c(NA,  varF.A)[interval.A]
        varF.B_times <- c(NA,  varF.B)[interval.B]

        ## compute variance of logFR
        varlogFR_times <- varF.A_times/(F.A_times^2) + 
                            varF.B_times/(F.B_times^2)

        ## CI for logFR_times (matrix form)
        logFR_times_loCI <- logFR_times - qnorm(1-alpha/2)*sqrt(varlogFR_times)
        logFR_times_upCI <- logFR_times + qnorm(1-alpha/2)*sqrt(varlogFR_times)
        

        ## I can't see how this next line of code could be necessary, but
        ## haven't had time to test thoroughly to ensure that it's not,
        ## so it remains for now...
        nas <- ( is.na(F.A_times) | is.na(F.B_times) )
        if ( any( nas ) )
            varlogFR_times[ nas ] <- NA

        ## create data.frame containing: point estimate, variance, 
        ## confidence interval, etc.
        outObj <-
            data.frame(
                evalTime  = times,
                nEvents   = NA,
                nEvents.1 = c(0, cumsum(trtA$n.event))[interval.A],
                nEvents.2 = c(0, cumsum(trtB$n.event))[interval.B],
                FR        = FR_times,  
                varlogFR  = varlogFR_times,
                FR_loCI   = exp( logFR_times_loCI ),
                FR_upCI   = exp( logFR_times_upCI )
              )

        ## fill in remaining components 
        outObj <- transform( outObj, 
                      nEvents = nEvents.1 + nEvents.2,
                      VE      = 1 - FR,   
                      VE_loCI = 1 - FR_upCI,
                      VE_upCI = 1 - FR_loCI )

        ## Check for 0 counts in groups. If either group has 0 infections
        ## then we will use 'VEci_binom' to compute a CI for the VE and
        ## fill it into 'outObj'.  
        if ( any( ind <- ( outObj$nEvents.1==0 | outObj$nEvents.2==0 ) )) {

            ci <- VEci_binom( 
                      N  = outObj$nEvents[ind],  
                      Nv = outObj$nEvents.1[ind],
                      randFraction = randFraction,
                      alpha = alpha )

            outObj$VE_loCI[ind] <- ci[ ,1]
            outObj$VE_upCI[ind] <- ci[ ,2]

            outObj$FR_loCI[ind] <- 1 - outObj$VE_upCI[ind]
            outObj$FR_upCI[ind] <- 1 - outObj$VE_loCI[ind]
        } 
        return( outObj )
    }

    force(times)

    ## Apply our 'cir' function to the list of datasets 'd', and return the list
    cirList <- lapply(d, FUN=cir, times=times, arms=arms, formula=survForm,
                      alpha=alphaLevel, randFraction = randFraction)

    ## If SIMPLIFY is turned off, return object and exit
    if ( !SIMPLIFY ) 
       return( cirList ) 

    ## If a single trial was passed via 'd' then cirList is a length 1 list
    ## containing a data.frame -  we remove it from the data.frame
    if (length(d) == 1) 
       return( cirList[[1]] )

    ## if a single time is provided, then rbind together all the one-row
    ##  data.frames into a single one
    if (length(times) == 1 && times != "all") 
       return( do.call(rbind, cirList) )

    ## if none of the criteria above were satisfied
    return( cirList ) 
}




## computes unadjusted hazard ratio and CI on given dataset(s)
##
## Not generic at all right now, assumes that:
##   - data comes from a std. trial data.frame
##   - that formula is: Surv(exit-entry, event) ~ trt
##
## Not sure how to generalize it as not sure where/how I'd use a general version

coxHR <- function(d, arms=c(1,0), alphaLevel=0.05, randFraction = 0.5,
                   ..., SIMPLIFY=TRUE)
{
    ##
    ##  d      a data.frame or list of data.frames to be used for fitting the model
    ##
    ##  arms   vector of length 2.   Indicates which treatment arms to use in
    ##         computing the hazard ratio.  The values should be the numbers (or
    ##         names, if assigned), for the treatments.  The first value of arms
    ##         specifies the treatment to be used in the numerator of the ratio,
    ##         and the 2nd value the treatment to be compared to (i.e. in the 
    ##         treatment used in the denominator).
    ##
    ## alphaLevel: the two-sided alpha-level used to construct confidence 
    ##             intervals for the estimate
    ##
    ## survForm: the formula to be used in survfit() to create summary info
    ##           used to compute the estimator.
    ##
    ## randFraction:  the fraction of randomization going to the particular
    ##                treatment group, when only that treatment group and the
    ##                placebo group are considered
    ##
    ##  ...    Not used at present.  Including this argument force exact 
    ##         specification of argument 'SIMPLIFY' if user desires to change
    ##         it (partial matching doesn't work after "..." arguments)
    ##
    ##  SIMPLIFY  Logical - Should the return structure be simplified from a 
    ##            list of data.frames to a single data.frame?
    ##
   
    survForm <- Surv(exit-entry, event) ~ trt
    
    ## force evaluation of objects so that they can be successfully passed
    ## through to other functions
    force(d)
    force(arms)
    force(randFraction)

    ## is 'd' is a single data.frame, convert it into a list 
    if (is.data.frame(d))  
      d <- list(d)


    ## restrict data to specified arms, if necessary 
    d <- lapply(d, function(ds,vals) ds[ ds$trt %in% vals, ], vals=arms )

    ## applies coxph to a single dataset
    cph <- function(d.i, formula, alpha, randFraction) {
        force(d.i)
        coxPH <- coxph(formula, data=d.i)

        ## extract the hazard ratio for 'trt' and store
        ## (note: all 'as.vector' calls are used to strip off names since they can
        ## cause the code to malfunction (e.g. in creation of data.frame - names
        ## override the names I'm trying to give the components!)
        logHR      <- as.vector( coef(coxPH) )
        varlogHR   <- as.vector( vcov(coxPH) )
        logHR_loCI <- logHR - qnorm(1 - alpha/2)*sqrt(varlogHR)
        logHR_upCI <- logHR + qnorm(1 - alpha/2)*sqrt(varlogHR)


        ## next code chunk is ancillary to the estimate, used only to
        ## retrieve info on number of infections in the groups
        SurvSumm <- summary( survfit(formula, data=d.i) )
        SurvSumm <- as.data.frame( SurvSumm[c("strata","time","n.event")] )

        ## if there are no events in either group, then punt.
        ## Fill in return object with NAs (mostly) and return
        if ( nrow(SurvSumm)==0 ) {
            return( data.frame(
                        nEvents   =  0,
                        nEvents.1 =  0,
                        nEvents.2 =  0,
                        HR        = NA,
                        varlogHR  = NA,
                        HR_loCI   = NA,
                        HR_upCI   = NA,
                        VE        = NA,
                        VE_loCI   = NA,
                        VE_upCI   = NA ) )
        }

        eventsByStrata <- split( SurvSumm, 
                                 factor( SurvSumm$strata,
                                         levels=paste0("trt=", arms ))) 

        events <-  sapply( eventsByStrata,
                           function(x) {
                             if (nrow(x)>0)
                               cumsum(x$n.event)[nrow(x)]
                             else  0 } )

        ## create data.frame containing: point estimate, variance, 
        ## confidence interval, etc.
        outObj <-
            data.frame(
                nEvents   = sum(events),
                nEvents.1 = events[1],
                nEvents.2 = events[2],
                HR        = exp(logHR),
                varlogHR  = varlogHR,
                HR_loCI   = exp( logHR_loCI ),
                HR_upCI   = exp( logHR_upCI ) 
            )

        ## fill in remaining components 
        outObj <- transform( outObj,
                      VE      = 1 - HR,
                      VE_loCI = 1 - HR_upCI,
                      VE_upCI = 1 - HR_loCI )

        ## Check for 0 counts in groups. If either group has 0 infections
        ## then we will use 'VEci_binom' to compute a CI for the VE and
        ## fill it into 'outObj'.  
        if ( any( ind <- ( outObj$nEvents.1==0 | outObj$nEvents.2==0 ) )) {

            ci <- VEci_binom(
                      N  = outObj$nEvents[ind],
                      Nv = outObj$nEvents.1[ind],
                      randFraction = randFraction,
                      alpha = alpha )

            outObj$VE_loCI[ind] <- ci[ ,1]
            outObj$VE_upCI[ind] <- ci[ ,2]

            outObj$HR_loCI[ind] <- 1 - outObj$VE_upCI[ind]
            outObj$HR_upCI[ind] <- 1 - outObj$VE_loCI[ind]
        }

        return( outObj )
    }

    ## Apply our 'cph' function to the list of datasets 'd', and return the list
    cphList <- lapply(d, FUN=cph, formula=survForm, alpha=alphaLevel,
                      randFraction = randFraction)

    ## If SIMPLIFY is turned off, return object and exit
    if ( !SIMPLIFY ) 
       return( cphList ) 

    ## If a single trial was passed via 'd' then cphList is a length 1 list
    ## containing a data.frame - we remove it from the list 
    if (length(d) == 1) 
       return( cphList[[1]] )

    ## else if multiple trial were provided, then rbind together the one-row
    ##  data.frames into a single one
    return( do.call(rbind, cphList) )
}



## This function implements a "final" analysis, just based on stage 1 data, 
## which will be used for two purposes:
##
##   (1) to perform an "End-of-Stage-1" test, when the final ppt. completes
##       their stage 1 follow-up (if that timepoint is reached without
##       encountering a stopping boundary first).  This test will determine
##       whether the trial continues to stage2, and also whether there is 
##       efficacy at the level of the design alternative (finalVE)
##
##   (2) to perform the "last-stage1-analysis" in cases when a stopping boundary
##       is hit before the end-of-stage 1 is reached.  Because this test is to
##       determine if there is evidence of efficacy, there is no need to do the
##       test if stopping was for harm or nonEfficacy.  So we do this test only
##       if we stop for highEff during stage 1.
finalStage1Test <- function(dat, 
                       analysisType=c("final","stopTime"), 
                       stage1VE=NULL,
                       lowerVE,
                       alphaLevel, 
                       estimand=c("cox","cuminc","combined"),
                       time=NULL, 
                       boundLabels=c("Eff", "NonEffFinal"),
                       randFraction = null.p )

## Arguments:
##   dat            a trial or list of trials
##
##   analysisType   Which type of analysis should be done
##
##   stage1VE       Used to determine whether the trial should move to stage 2.
##                  Only applicable if analysisType == "Final".
##                  The lower bound of the VE confidence interval must exceed
##                  this value.  (Minimal VE threshhold for further study)
##
##   lowerVE        vector of VE values to test (using stage 1 data only) when
##                  the trial stops, either by hitting a bound of by reaching 
##                  the end of stage 1.
##
##   alphaLevel     alpha level to use in constructing the confidence interval
##                  for the VE estimate that will be compared to 'lowerVE' and
##                  (if applicable) 'stage1VE'.
##
##   estimand       which estimand(s) to use to estimate VE
##
##   time           (Only used if analysisType == "stopTime") 
##                  The calendar time at which to censor the stage1 data prior
##                  to estimating VE and comparing to 'finalVE'
##
##   boundLabels    (Used only if analysisType == "final").  
##                  These are values to return for the result component 
##                  "boundHit", when the comparison to 'lowerVE' is TRUE,
##                  or FALSE respectively.  
##                  Call the lower bound for the CI of the VE estimates 'loCI':
##                     if loCI > lowerVE for *all* estimands specified,                
##                          then boundHit <- boundLabels[1]
##                     if loCI <= lowerVE, for *any* estimand specified
##                          then  boundHit <-  boundLabels[2]
##
##  randFraction    Randomization fraction for current treatment group relative
##                  to that of the placebo group.  This information is used only
##                  in the case of having 0 infections in one of the groups
##                  (active or placebo), so that a proper confidence interval 
##                  cannot be constructed.  In that case, we pretend that our
##                  VE estimate is:  1 - (Nvacc/Nplac)*(Rp/Rv)
##                  and that Nvacc ~ binomial(Nvacc+Nplac, p).
##                  We compute an exact CI for 'p' based on the Clopper-Pearson
##                  method and use it to compute a CI for VE as just defined.
{
    force( randFraction )

    analysisType  <- match.arg(analysisType)
    estimand      <- match.arg(estimand)

    ## Set indicator of end of stage 1 analysis
    EndStg1 <- (analysisType == "final")

    ## force evaluation of a few arguments
    force(alphaLevel)
    force(lowerVE)
    if ( !missing(stage1VE) ) 
        force(stage1VE)
    if ( !missing(time) ) 
        force(time)

    ## TEMP CODE - need to modify monitorTrials to specify a specific estimand
    ## (not "combined") to use for the final stage 1/stop-time  analysis
    if (estimand == "combined")
       estimand <- "cuminc"

    ## set indicators for estimation types to use
    cox  <- (estimand %in% c("cox",   "combined"))
    cir  <- (estimand %in% c("cuminc","combined"))

    if (analysisType=="stopTime") {
        if (is.null(time)) {
            stop("Argument 'time' must be specified when analysisType='stopTime'\n")
        }
        ## censor data to 'time'
        dat <- censorTrial(dat, times=time, timeScale="calendar")
    }


    if (cir) {
        cumIncOut <- cumIncRatio( dat, 
                         times="last", 
                         alphaLevel = alphaLevel,
                         randFraction = randFraction )

        CIRobj <- transform( cumIncOut,
                      infectTotal = nEvents,
                      infectSplit = paste0("Pl:Vx =", nEvents.2, ":", nEvents.1),
                      nPlac       = nEvents.2,
                      nVacc       = nEvents.1,
                      VE_CIR      = VE,
                      VE_CIR_loCI = VE_loCI,
                      VE_CIR_upCI = VE_upCI )

        if (EndStg1) {
            ## set efficacy result based on lower bound of CI, unless we have no
            ## CI because we have no infections in vacc. group, then set to TRUE
            effResult <- ifelse(CIRobj$VE_loCI > stage1VE, TRUE, FALSE )
        }

        altDetected_CIR <- ifelse(CIRobj$VE_loCI > lowerVE, TRUE, FALSE)

             
        ## subset to retain only the var.s created above
        CIRobj <- CIRobj[,c("infectTotal", "infectSplit", "nPlac", "nVacc",
                           "VE_CIR", "VE_CIR_loCI", "VE_CIR_upCI") ]
    }

    if (cox) {
        coxHRout <- coxHR( dat, 
                           alphaLevel = alphaLevel,
                           randFraction = randFraction )

        CoxObj <- transform( coxHRout,
                      infectTotal = nEvents,
                      infectSplit = paste0("Pl:Vx =", nEvents.2, ":", nEvents.1),
                      nPlac       = nEvents.2,
                      nVacc       = nEvents.1,
                      VE_Cox      = VE,
                      VE_Cox_loCI = VE_loCI,
                      VE_Cox_upCI = VE_upCI)

        if (EndStg1) {
            ## set efficacy result based on lower bound of CI, unless we have no
            ## CI because we have no infections in vacc. group, then set to TRUE
            effResult <- ifelse(CoxObj$VE_loCI > stage1VE, TRUE, FALSE )
        }

        altDetected_Cox <- ifelse( CoxObj$VE_loCI > lowerVE, TRUE, FALSE )


        ## subset to retain only the var.s created above
        CoxObj <- CoxObj[,c("infectTotal", "infectSplit", "nPlac", "nVacc",
                           "VE_Cox",     "VE_Cox_loCI", "VE_Cox_upCI") ]
    }

    if (cox && cir) {
        altDetected <- altDetected_CIR & altDetected_Cox
        ## cbind the output together, dropping duplicated columns
        summObj <- cbind( CIRobj, 
                          CoxObj[, !names(CoxObj) %in% names(CIRobj) ] ) 
    } else {
        if (cox) {
            altDetected <- altDetected_Cox
            summObj     <- CoxObj
        } else {
            altDetected <- altDetected_CIR
            summObj     <- CIRobj
        }
    } 

    names( altDetected ) <- as.character( lowerVE )

    if ( EndStg1 ) {
        
        ## using 'as.vector' to strip the names off of 'altDetected' so they don't
        ## transfer onto 'boundHit;
        boundHit <- ifelse(effResult, boundLabels[1], boundLabels[2])

        outList <- 
            list( boundWasHit = TRUE,
                  boundHit = boundHit,
                  boundType = "FinalStage1",
                  ## only assign stopTime for nonefficacy result
                  stopTime = ifelse(effResult, NA, max( dat$exit)),
                  stage1summ = summObj, 
                  altDetected = altDetected )
    } else {
        ## this analysis comes after a bound was already hit, so we're just 
        ## reporting out on the last analysis.  So we won't output things
        ## like 'boundHit', and stopTime, as those have already been defined
        ## We just output the summaryObject with estimates and CIs and 
        outList <- list( stage1summ = summObj,
                         altDetected = altDetected )
    }
    return( outList )
}


## This function is used to create a confidence interval for the VE when
## we have 0 infections in either the active or placebo groups. It is
## *ONLY* for that case - no other.  It is used because we cannot construct
## a 'usual' CI in this setting.
##
## The approach is based off estimating VE as:   
##   VE = 1 - (Nv/Np)*(Rp/Rv)  
##      with Nv, Np being infection counts in the vacc. and plac. groups
##      and Rp, Rv the randomization fractions of those grps. 
##
## We treat Nvacc as ~ binomial(N, p), where N = (Nv+Np),  and compute an
## exact confidence interval for p using the Clopper-Pearson approach.
## We then use that CI to derive a CI for VE, based on the equation given
## a few lines above here.
##
##  N         = total infections 
##  Nv        = number of infections in vaccinee group
##  alpha     = the 2-sided confidence level
##  randFraction = randomization fraction of the current active trt group relative
##              to active + placebo.
## 
VEci_binom <- function(N, Nv, randFraction=0.5, alpha=0.05) {

    if ( any( (Nv != 0  &  Nv != N) ) )
        stop("This function is only usable when Nv equals 0 or N \n\n")
    
    L <- length( N )
    lowerBound <- numeric( L )
    upperBound <- numeric( L )

    ## logical indicators of 0 and N infections
    ind0 <- (Nv == 0) 
    indN <- (Nv == N)

    ## Bounds first for binomial 'p'
    lowerBound[ ind0 ] <- 0
    upperBound[ ind0 ] <- 1 - (alpha/2)^(1/N)

    lowerBound[ indN ] <- (alpha/2)^(1/N)
    upperBound[ indN ] <- 1 
 
    ## This is the reversed CI for binomial 'p', the make
    ## next set of computations simpler.
    revCI <- cbind(upperBound, lowerBound)

    ## Then transform it to get CI for VE
    ## First re-write VE equation as:  
    ##   VE = 1 - (Nv/N)/((N-Nv)/N)*Rp/Rv =   1 - [p/(1-p)] * Rp/Rv
    ## 
    ## then let (pLL, pUL) = lower-limit and upper-limit of the CI for 'p'
    ## and note:
    ##  since - p/(1-p) is decreasing in 'p', the lower/upper CI limits 
    ##    for the CI of -p/(1-p) will be: (-pUL/(1-pUL), -pLL/(1-pLL))  
    ##    and the CI for VE comes by direct manipulation of this
   
    ci <- cbind(lowerBound, upperBound)

    VE_ci <- 1 - revCI/(1-revCI) * (1-randFraction)/randFraction

    return( VE_ci )
}


## Function that, given a vector of infection counts and a simulated trial 
## dataset 'd' will return the calendar times at which those infection counts
## were reached.  The times are needed for censoring the data prior to 
## computing various statistics for, e.g. for performing non-efficacy or 
## high-efficacy analyses.


#* Suggestion to think about:
#*    can we unlink 'cnts' from length of d?  Then we'd be able to 
#*    allow for 'cnts' to be unspecified so that we can get the times 
#*    for every infection in d (or d[[j]])
#*
#*    Not sure what bad things can happen if we allow that...
#*    or if it makes one of the other usages confusing... 


## Arguments:
##
##    d  either a trial data.frame (i.e. data from a single trial) or a list
##       of such data.frames (each representing a trial)
##
## cnts: integer vector giving the infection counts (summed over the treatment
##       arms given by argument 'arms') for which associated calendar times
##       are desired. 
## 
##       If 'cnts' is a vector, the values will be applied to every trial
##
##       If 'cnts' is a list, then each component should contain a vector 
##       of counts for a different trial.  i.e. cnts[[j]] will contain the
##       counts for which we want times for trial d[[j]].  Because of this
##       relationship between 'd' and 'cnts', both must be the same length.
##
## arms:  The treatment arm codes for the ppt.s whose infections should be 
##        counted when determining the infection totals.  If unspecified
##        all data is used (i.e. arms are not considered)
##
##  ... :  Not used.  Purpose is to force exact specification of arguments
##         that follow it.
##
## SIMPLIFY:  Logical value indicating whether the function should return 
##            a vector rather than a length 1 list containing the vector
##            when 'd' is a single data.frame (or list of length 1)
##
##
## Return Value:
##
##   A list the same length as 'd' (length=1 if 'd' is a data.frame), each
##   component of which contains a vector of calendar times at which the 
##   associated infection totals (given by, or in 'cnts') is reached.
##

##getInfectionTimes <- function( d , cnts, arms=c(0,1), ..., SIMPLIFY=TRUE)
getInfectionTimes <- function( d , cnts, arms=NULL, ..., SIMPLIFY=TRUE)
{

  ## if 'cnts' is a list, then 'd' must be too
  if (is.list(cnts) && (!is.list(d) || length(d)!=length(cnts)) )
      stop("When 'cnts' is a list, then 'd' must also be a list and their",
           " lengths\n", "must be the same\n\n")

  ## coerce 'cnts' to be integer valued (within the list, if present)
  cnts <- if (is.list(cnts)) {
              lapply(cnts, as.integer)
          } else {
              as.integer(cnts)
          } 
  ## place 'd' inside a list if it's just a single data.frame
  if ( is.data.frame(d) )
     d <- list( d )

  ## define function to 'apply' to all datasets (lapply or mapply)
  laf <- function(x, cnt, arms=NULL) {

      cnt.int <- as.integer(cnt)
      if ( is.null(arms) ) {
          eventTimes <- sort( x$exit[ x$event == 1 ])
      } else {
          eventTimes <- sort( x$exit[ x$trt %in% arms & x$event == 1 ])
      }

      ## restrict to counts that were attained in this trial
      cnt.obt <- cnt.int[ cnt.int <= length(eventTimes) ]

      ## return times corresponding to the count totals obtained
      eventTimes[ cnt.obt ]  
  }

  Times <- if (is.list(cnts)) {
              mapply(FUN=laf, x=d, cnt=cnts)
          } else {
              lapply(X=d, FUN=laf, cnt=cnts)
          }

  if (length(d) > 1 || !SIMPLIFY) { 
      return( Times )
  } else {
      return( Times[[1]] )
  }
}



## Censors trial data.frames to specified times:
##
##   d:  either a trial data.frame (i.e. data from a single trial) or a list
##       of such data.frames (each representing a trial)
##
##  times: gives the time(s) to use for censoring.  
## 
##         If 'times' is a vector, then the following behavior will be obtained:
##            every trial will be censored using every value of 'times'.
##
##         If 'times' is a list, then it must have length equal to length(d),
##         and the resultant behavior will be to apply the times from times[[j]]
##         to trial d[[j]] - i.e. each component of 'times' specifies a set of
##         times to use in censoring a particular trial.
## 
##  arms:  UNDER CONSTRUCTION - NOT CURRENTLY IN USE - WILL BE IGNORED.
##         vector of codes of the treatment arms that should be censored.
##         Default is to censor all arms.
##
##  timeScale:  Specifies whether the censoring should be done on calendar time 
##         or follow-up time - i.e. whether the values in 'times' should be taken
##         as time since the study began (calendar time), or the time since each
##         participant enrolled (follow-up time).
##
##  type:  Type of censoring to perform. Possible values are 'right' and 'left'
##         The default is 'right'.  Right censoring discards data accrued after
##         a given time, whereas left censoring discards data accrued *before*
##         a given time.  Interval censoring can be achieved by employing the
##         function twice, once with right censoring and once with left. 
##
##   ...:  Not currently used
##
##  SIMPLIFY:  Indicates whether a data.frame should be returned instead of a 
##             list of length 1, in the case where a single data.frame (or list 
##             length 1) and single time are provided by the user.
##             Or, in the case of mulitple trials and a single time per trial,
##             whether the nested list structure should be simplified to a 
##             un-nested list (i.e. the inner list be "unlist()ed")

censorTrial <- function(d, times, arms=NULL, timeScale=c("calendar","follow-up"),
                        type=c("right","left"), ..., SIMPLIFY=TRUE)
{
    ## To handle the case when 'd' is a list of trials, we force it to always
    ## be a list (i.e. we coerce to list when needed)
    if ( is.data.frame(d) )
        d <- list(d)

    ## determine how many trials we're censoring
    nTrials <- length(d)
    nTimes  <- length(times) 

    if (is.list(times) && nTimes != nTrials)
      stop("Argument 'times' is a list and has length different than 'd' - ",
           "this is not permitted\n")

    ## if 'arms' has been specified, subset data using it
    #if (!is.null(arms)) {
    #  d <- lapply(d, function(x) x[x$trt %in% arms, ] )
    #}
             
    timeScale <- match.arg(timeScale)
    type <- match.arg(type)


    ### Define 4 censoring functions, for all combinations of type and timeScale

    ## (1) function to right censor a single trial at a single calendar time
    right.censor.Calendar <- function(trial, time) {
          ## restrict to those enrolled before 'time'
          cens <- trial[trial$entry<time, ]

          ## censor events that occur after 'time'
          cens$event[ cens$exit > time ] <- 0

          ## alter 'exit' time to reflect censoring
          cens$exit <- pmin( cens$exit, time)

          return(cens)
    }

    ## (2) function to left censor a single trial at a single calendar time
    left.censor.Calendar <- function(trial, time) {
          ## restrict to those still active at time 'time'
          cens <- trial[trial$exit > time, ]

          ## alter 'entry' time to reflect censoring, it must be >= time
          cens$entry <- pmax( cens$entry, time)

          return(cens)
    }

    ## (3) function to right censor a single trial at a single follow-up time
    ## follow-up time = (exit_time - entry_time)
    right.censor.Followup <- function(trial, time) {

          ## No exclusion of participants
          cens <- trial

          ## censor events with follow-up time greater than 'time'
          cens$event[ (cens$exit - cens$entry) > time ] <- 0

          ## alter 'exit' time to reflect censoring
          cens$exit <- pmin( cens$exit, cens$entry + time ) 

          return(cens)
    }

    ## (4) function to left censor a single trial at a single follow-up time
    left.censor.Followup <- function(trial, time) {
          ## restrict to participants with at least 'time' followup
          cens <- trial[ (trial$exit - trial$entry) > time, ]

          ## -- not sure that we should do this -- leave it out for now...
          ## alter 'entry' time to reflect censoring (move up by 'time') to
          ## reduce FU by 'time' 
          ## cens$entry <- cens$entry + time 

          return(cens)
    }
    
    ## choose appropriate censoring function
    censFunc <- switch( type, 
                  "right"= switch(timeScale, "calendar"=right.censor.Calendar,
                                    "follow-up"=right.censor.Followup, NA ),
                  "left" = switch(timeScale, "calendar"=left.censor.Calendar,
                                    "follow-up"=left.censor.Followup, NA ) )


    ## apply censoring
    if (!is.list(times)) {
        ## apply all times to all trials.  This returns a nested list:
        ##   outer list over trials, inner list over times within trial
        cens <- lapply(d, FUN = function(d.i, time, cenFunc) 
                         lapply(X=time, FUN=cenFunc, trial=d.i),
                       time=times, cenFunc=censFunc )
    } else {
        ## otherwise if 'times' is a list (same length as 'd') do: 
        cens <- mapply(FUN= function(d.i, time, cenFunc) 
                          lapply(X=time, FUN=cenFunc, trial=d.i ),
                       d.i=d, time=times, 
                       MoreArgs=list(cenFunc=censFunc),
                       SIMPLIFY=FALSE, USE.NAMES=FALSE)
    }

  if ( SIMPLIFY ) {
    ## - If the output is only one data.frame, remove it from its list,
    ## - If the output is one trial at mulitple times, remove the outer list.
    ## - If the output is multiple trials each at a single time, then remove
    ##     the inner list.
    ## [note: 'unlist()' not used as it strips off the data.frame attribute ]

    if (nTrials==1 && nTimes==1)
        return( cens[[1]][[1]] )

    if (nTrials==1 && nTimes>1)
        return( cens[[1]] )

    if (nTrials>1 && (nTimes==1 || 
          (is.list(times) && all(unlist(lapply(times,length))==1)) ) )
        return( lapply(cens, function(x) x[[1]] ) )

    return( cens )

  } else {
    return ( cens )
  }
}



#################   Begin code for monitorTrial  ################

## Some arguments described just below the argument list 

#' Group Sequential Monitoring of Simulated Efficacy Trials for the Event of Potential Harm, Non-Efficacy, and High Efficacy
#'
#' \code{monitorTrial} applies a group sequential monitoring procedure to data-sets generated by \code{simTrial}, which may result in modification or termination of each simulated trial.
#'
#' @param dataFile if \code{saveDir = NULL}, a list returned by \code{simTrial}; otherwise a name (character string) of an \code{.RData} file created by \code{simTrial}
#' @param stage1 the final week of stage 1 in a two-stage trial
#' @param stage2 the final week of stage 2 in a two-stage trial, i.e., the maximum follow-up time
#' @param harmMonitorRange a 2-component numeric vector specifying the range of the pooled number of infections (pooled over the placebo and vaccine arm accruing infections the fastest) over which the type I error rate, specified in \code{harmMonitorAlpha}, will be spent (per vaccine arm). Note that \code{harmMonitorRange} does not specify a range for which potential-harm stopping boundaries will be computed; instead, it specifies when potential-harm monitoring will start, and the range over which \code{harmMonitorAlpha} will be spent.
#' @param harmMonitorAlpha a numeric value (0.05 by default) specifying the overall type I error rate for potential-harm monitoring (per vaccine arm). To turn off potential-harm monitoring, set \code{harmMonitorAlpha} equal to 0.00001.
#' @param alphaPerTest a per-test nominal/unadjusted alpha level for potential-harm monitoring. If \code{NULL}, a per-test alpha level is calculated that yields a cumulative alpha of \code{harmMonitorAlpha} at the end of \code{harmMonitorRange}.
#' @param nonEffStartMethod a character string specifying the method used for determining when non-efficacy monitoring is to start. The default method of Freidlin, Korn, and Gray (2010) ("\code{FKG}") calculates the minimal pooled infection count (pooled over the placebo and vaccine arm accruing infections the fastest) such that a hazard-ratio-based VE point estimate of 0\% would result in declaring non-efficacy, i.e., the upper bound of the two-sided (1-\code{alphaNoneff}) x 100\% confidence interval for VE based on the asymptotic variance of the log-rank statistic equals the non-efficacy threshold specified as component \code{upperVEnonEff} in the list \code{nonEffStartParams}. If this list component is left unspecified, the argument \code{upperVEnonEff} is used as the non-efficacy threshold. The alternative method ("\code{fixed}") starts non-efficacy monitoring at a fixed pooled infection count (pooled over the placebo and vaccine arm accruing infections the fastest) specified by component \code{N1} in the list \code{nonEffStartParams}.
#' @param nonEffStartParams a list with named components specifying parameters required by \code{nonEffStartMethod} (\code{NULL} by default)
#' @param nonEffIntervalUnit a character string specifying whether intervals between two adjacent non-efficacy interim analyses should be event-driven (default option "\code{counts}") or calendar time-driven (option "\code{time}")
#' @param nonEffInterval a numeric vector (a number of infections or a number of weeks) specifying the timing of non-efficacy interim analyses. If a single numeric value is specified, then all interim looks are equidistant (in terms of the number of infections or weeks), and the value specifies the constant increment of information or time between two adjacent interim looks. If a numeric vector with at least two components is specified, then, following the initial interim look, the timing of subsequent interim looks is determined by (potentially differential) increments of information or time specified by this vector.
#' @param lowerVEnoneff specifies criterion 1 for declaring non-efficacy: the lower bound of the two-sided (1-\code{alphaNoneff}) x 100\% confidence interval(s) for the VE estimand(s) lie(s) below \code{lowerVEnoneff} (typically set equal to 0). If \code{NULL} (default), this criterion is ignored.
#' @param upperVEnoneff specifies criterion 2 for declaring non-efficacy: the upper bound of the two-sided (1-\code{alphaNoneff}) x 100\% confidence interval(s) for the VE estimand(s) lie(s) below \code{upperVEnoneff} (typically a number in the 0--0.5 range)
#' @param highVE specifies a criterion for declaring high-efficacy: the lower bound of the two-sided (1-\code{alphaHigh}) x 100\% confidence interval for the VE estimand lies above \code{highVE} (typically a number in the 0.5--1 range). To turn off high efficacy monitoring, set \code{highVE} equal to 1.
#' @param stage1VE specifies a criterion for advancement of a treatment's evaluation into Stage 2: the lower bound of the two-sided (1-\code{alphaStage1}) x 100\% confidence interval for the VE estimand lies above \code{stage1VE} (typically set equal to 0)
#' @param lowerVEuncPower a numeric vector with each component specifying a one-sided null hypothesis H0: VE(0--\code{stage1}) \eqn{\le} \code{lowerVEuncPower} x 100\%. Unconditional power (i.e., accounting for sequential monitoring) to reject each H0 is calculated, where the rejection region is defined by the lower bound of the two-sided (1-\code{alphaUncPower}) x 100\% confidence interval for the VE estimand being above the respective component of \code{lowerVEuncPower} (typically values in the 0--0.5 range).
#' @param alphaNoneff one minus the nominal confidence level of the two-sided confidence interval used for non-efficacy monitoring
#' @param alphaHigh one minus the nominal confidence level of the two-sided confidence interval used for high efficacy monitoring
#' @param alphaStage1 one minus the nominal confidence level of the two-sided confidence interval used for determining whether a treatment's evaluation advances into Stage 2
#' @param alphaUncPower one minus the nominal confidence level of the two-sided confidence interval used to test one-sided null hypotheses H0: VE(0-\code{stage1}) \eqn{\le} \code{lowerVEuncPower} x 100\% against alternative hypotheses H1: VE(0--\code{stage1}) \eqn{>} \code{lowerVEuncPower} x 100\%. The same nominal confidence level is applied for each component of \code{lowerVEuncPower}.
#' @param estimand a character string specifying the choice of VE estimand(s) used in non- and high efficacy monitoring, advancement rule for Stage 2, and unconditional power calculations. Three options are implemented: (1) the `pure' Cox approach (\code{"cox"}), where VE is defined as 1-hazard ratio (treatment/control) and estimated by the maximum partial likelihood estimator in the Cox model; (2) the `pure' cumulative incidence-based approach (\code{"cuminc"}), where VE is defined as 1-cumulative incidence ratio (treatment/control) and estimated by the transformation of the Nelson-Aalen estimator for the cumulative hazard function; and (3) the combined approach (\code{"combined"}), where both aforementioned VE estimands are used for non-efficacy monitoring while the cumulative VE estimand is used for all other purposes. Only the first three characters are necessary.
#' @param laggedMonitoring a logical value (\code{FALSE} by default) indicating whether "per-protocol" non-efficacy monitoring should additionally be conducted for events occurring after \code{lagTime} weeks as a more conservative non-efficacy monitoring approach. If \code{TRUE} and \code{estimand = "combined"}, the cumulative VE estimand is considered only for non-efficacy monitoring.
#' @param lagTime a time point (in weeks) defining the per-protocol VE estimand, i.e., VE(\code{lagTime}--\code{stage1}). This VE estimand is also used in "per-protocol" non-efficacy monitoring if \code{laggedMonitoring} equals \code{TRUE}. It is typically chosen as the date of the last immunization or the date of the visit following the last immunization.
#' @param saveFile a character string specifying the name of the output \code{.RData} file. If \code{NULL} (default), a default file name will be used.
#' @param saveDir a character string specifying a path for \code{dataFile}. If supplied, the output is also saved as an \code{.RData} file in this directory; otherwise the output is returned as a list.
#' @param verbose a logical value indicating whether information on the output directory, file name, and monitoring outcomes should be printed out (default is \code{TRUE})
#'
#' @details All time variables use week as the unit of time. Month is defined as 52/12 weeks.
#'
#' Potential harm monitoring starts at the \code{harmMonitorRange[1]}-th infection pooled over the placebo group and the vaccine regimen that accrues infections the fastest. The potential harm analyses continue at each additional infection until the first interim analysis for non-efficacy. The monitoring is implemented with exact one-sided binomial tests of H0: \eqn{p \le p0} versus H1: \eqn{p > p0}, where \eqn{p} is the probability that an infected participant was assigned to the vaccine group, and \eqn{p0} is a fixed constant that represents the null hypothesis that an infection is equally likely to be assigned vaccine or placebo. Each test is performed at the same prespecified nominal/unadjusted alpha-level (\code{alphaPerTest}), chosen based on simulations such that, for each vaccine regimen, the overall type I error rate by the \code{harmMonitorRange[2]}-th arm-pooled infection (i.e., the probability that the potential harm boundary is reached when the vaccine is actually safe, \eqn{p = p0}) equals \code{harmMonitorAlpha}.
#' 
#' Non-efficacy is defined as evidence that it is highly unlikely that the vaccine has a beneficial effect measured as VE(0--\code{stage1}) of \code{upperVEnoneff} x 100\% or more. The non-efficacy analyses for each vaccine regimen will start at the first infection (pooled over the vaccine and placebo arm) determined by \code{nonEffStartMethod}. Stopping for non-efficacy will lead to a reported two-sided (1-\code{alphaNoneff}) x 100\% CI for VE(0--\code{stage1}) with, optionally, the lower confidence bound below \code{lowerVEnoneff} and the upper confidence bound below \code{upperVEnoneff}, where \code{estimand} determines the choice of the VE(0--\code{stage1}) estimand. This approach is similar to the inefficacy monitoring approach of Freidlin, Korn, and Gray (2010). If \code{estimand = "combined"}, stopping for non-efficacy will lead to reported (1-\code{alphaNoneff}) x 100\% CIs for both VE parameters with, optionally, lower confidence bounds below \code{lowerVEnoneff} and upper confidence bounds below \code{upperVEnoneff}. If \code{laggedMonitoring = TRUE}, stopping for non-efficacy will lead to reported (1-\code{alphaNoneff}) x 100\% CIs for both VE(0--\code{stage1}) and VE(\code{lagTime}--\code{stage1}) with, optionally, lower confidence bounds below \code{lowerVEnoneff} and upper confidence bounds below \code{upperVEnoneff}.
#' 
#' High efficacy monitoring allows early detection of a highly protective vaccine if there is evidence that VE(0--\code{stage2}) \eqn{>} \code{highVE} x 100\%. It is synchronized with non-efficacy monitoring during Stage 1, and a single high-efficacy interim analysis during Stage 2 is conducted halfway between the end of Stage 1 and the end of the trial. While monitoring for potential harm and non-efficacy restricts to \code{stage1} infections, monitoring for high efficacy counts all infections during \code{stage1} or \code{stage2}, given that early stopping for high efficacy would only be warranted under evidence for durability of the efficacy.
#' 
#' The following principles and rules are applied in the monitoring procedure:
#' \itemize{
#'   \item Exclude all follow-up data from the analysis post-unblinding (and include all data pre-unblinding).
#'   \item The monitoring is based on modified ITT analysis, i.e., all subjects documented to be free of the study endpoint at baseline are included and analyzed according to the treatment assigned by randomization, ignoring how many vaccinations they received (only pre-unblinding follow-up included).
#'   \item If a vaccine hits the harm boundary, immediately discontinue vaccinations and accrual into this vaccine arm, and unblind this vaccine arm (continue post-unblinded follow-up until the end of Stage 1 for this vaccine arm).  
#'   \item If a vaccine hits the non-efficacy boundary, immediately discontinue vaccinations and accrual into this vaccine arm, keep blinded and continue follow-up until the end of Stage 1 for this vaccine arm. 
#'   \item If and when the last vaccine arm hits the non-efficacy (or harm) boundary, discontinue vaccinations and accrual into this vaccine arm, and unblind (the trial is over, completed in Stage 1).
#'   \item Stage 1 for the whole trial is over on the earliest date of the two events: (1) all vaccine arms have hit the harm or non-efficacy boundary; and (2) the last enrolled subject in the trial reaches the final \code{stage1} visit.
#'   \item Continue blinded follow-up until the end of Stage 2 for each vaccine arm that reaches the end of \code{stage1} with a positive efficacy (as defined by \code{stage1VE}) or high efficacy (as defined by \code{highVE}) result.
#'   \item If at least one vaccine arm reaches the end of \code{stage1} with a positive efficacy or high efficacy result, continue blinded follow-up in the placebo arm until the end of Stage 2.
#'   \item Stage 2 for the whole trial is over on the earliest date of the two events: (1) all subjects in the placebo arm and each vaccine arm that registered efficacy or high efficacy in \code{stage1} have failed or been censored; and (2) all subjects in the placebo arm and each vaccine arm that registered efficacy or high efficacy in \code{stage1} have completed the final \code{stage2} visit.
#' }
#' 
#' The above rules have the following implications:
#' \itemize{
#'   \item If a vaccine hits the non-efficacy boundary but Stage 1 for the whole trial is not over, then one includes in the analysis all follow-up through the final \code{stage1} visit for that vaccine regimen, including all individuals accrued up through the date of hitting the non-efficacy boundary (which will be the total number accrued to this vaccine arm).
#'   \item If a vaccine hits the harm boundary, all follow-up information through the date of hitting the harm boundary is included for this vaccine; no follow-up data are included after this date.
#'   \item If and when the last vaccine arm hits the non-efficacy (or harm) boundary, all follow-up information through the date of hitting the non-efficacy (or harm) boundary is included for this vaccine; no follow-up data are included after this date.
#' }
#' 
#' @return If \code{saveDir} (and, optionally \code{saveFile}) is specified, the output list (named \code{out}) is saved as an \code{.RData} file in \code{saveDir} (the path to \code{saveDir} is printed); otherwise it is returned. The output object is a list of length equal to the number of simulated trials, each of which is a list of length equal to the number of treatment arms, each of which is a list with (at least) the following components:
#' \itemize{
#'   \item \code{boundHit}: a character string stating the monitoring outcome in this treatment arm, i.e., one of \code{"Harm"}, \code{"NonEffInterim"}, \code{"NonEffFinal"}, \code{"Eff"}, or \code{"HighEff"}. The first four outcomes can occur in Stage 1, whereas the last outcome can combine data over Stage 1 and Stage 2.
#'   \item \code{stopTime}: the time of hitting a stopping boundary since the first subject enrolled in the trial
#'   \item \code{stopInfectCnt}: the pooled number of infections at \code{stopTime}
#'   \item \code{summObj}: a \code{data.frame} containing summary information from each non-/high efficacy interim analysis
#'   \item \code{finalHRci}: the final CI for the hazard ratio, available if \code{estimand!="cuminc"} and there is at least 1 infection in each arm
#'   \item \code{firstNonEffCnt}: the number of infections that triggered non-efficacy monitoring (if available)
#'   \item \code{totInfecCnt}: the total number of \code{stage1} (\code{stage2} if \code{boundHit = "HighEff"}) infections
#'   \item \code{totInfecSplit}: a table with the numbers of \code{stage1} (\code{stage2} if \code{boundHit = "HighEff"}) infections in the treatment and control arm
#'   \item \code{lastExitTime}: the time between the first subject's enrollment and the last subject's exiting from the trial
#' }
#' 
#' @references Freidlin B., Korn E. L., and Gray R. (2010), A general inefficacy interim monitoring rule for randomized clinical trials. \emph{Clinical Trials} 7(3):197-208.
#'
#' @examples 
#' simData <- simTrial(N=c(1000, rep(700, 2)), aveVE=seq(0, 0.4, by=0.2), 
#'                     VEmodel="half", vePeriods=c(1, 27, 79), enrollPeriod=78, 
#'                     enrollPartial=13, enrollPartialRelRate=0.5, dropoutRate=0.05, 
#'                     infecRate=0.04, fuTime=156, 
#'                     visitSchedule=c(0, (13/3)*(1:4), seq(13*6/3, 156, by=13*2/3)),
#'                     missVaccProb=c(0,0.05,0.1,0.15), VEcutoffWeek=26, nTrials=5, 
#'                     stage1=78, randomSeed=300)
#'    
#' monitorData <- monitorTrial(dataFile=simData, stage1=78, stage2=156, 
#'                             harmMonitorRange=c(10,100), alphaPerTest=NULL, 
#'                             nonEffStartMethod="FKG", nonEffInterval=20, 
#'                             lowerVEnoneff=0, upperVEnoneff=0.4, highVE=0.7, 
#'                             stage1VE=0, lowerVEuncPower=0, alphaNoneff=0.05,
#'                             alphaHigh=0.05, alphaStage1=0.05, alphaUncPower=0.05,
#'                             estimand="cuminc", lagTime=26)
#'    
#' ### alternatively, to save the .RData output file (no '<-' needed):
#' ###
#' ### simTrial(N=c(1400, rep(1000, 2)), aveVE=seq(0, 0.4, by=0.2), VEmodel="half", 
#' ###          vePeriods=c(1, 27, 79), enrollPeriod=78, enrollPartial=13, 
#' ###          enrollPartialRelRate=0.5, dropoutRate=0.05, infecRate=0.04, fuTime=156, 
#' ###          visitSchedule=c(0, (13/3)*(1:4), seq(13*6/3, 156, by=13*2/3)), 
#' ###          missVaccProb=c(0,0.05,0.1,0.15), VEcutoffWeek=26, nTrials=30, 
#' ###          stage1=78, saveDir="./", randomSeed=300)
#' ###
#' ### monitorTrial(dataFile=
#' ###              "simTrial_nPlac=1400_nVacc=1000_1000_aveVE=0.2_0.4_infRate=0.04.RData", 
#' ###              stage1=78, stage2=156, harmMonitorRange=c(10,100), alphaPerTest=NULL, 
#' ###              nonEffStartMethod="FKG", nonEffInterval=20, lowerVEnoneff=0, 
#' ###              upperVEnoneff=0.4, highVE=0.7, stage1VE=0, lowerVEuncPower=0, 
#' ###              alphaNoneff=0.05, alphaHigh=0.05, alphaStage1=0.05, alphaUncPower=0.05, 
#' ###              estimand="cuminc", lagTime=26, saveDir="./")
#'
#' @seealso \code{\link{simTrial}}, \code{\link{censTrial}}, and \code{\link{rankTrial}}
#'
#' @export
monitorTrial <- function (dataFile,
                          stage1,
                          stage2,
                  
                          ## range over which to "spend" the type-I error specified in argument
                          ## 'harmMonitorRange'.  It should be a vector of length 2 giving the 
                          ## (start, stop) infection-count range over which the type-I error 
                          ## will be spent.  Please note that this argument does *NOT* dictate
                          ## when harm monitoring will END.  The range dictates when it will
                          ## START, and over what range the type-I error will be spread.  If
                          ## you plan to use nonEffStartMethod "FKG" or "fixed" then the 
                          ## 'stop' value of the range, if provided by the user, will not be
                          ## used and need not be specified. It may also be replaced with an NA.
                  
                          ## bounds are created.  The 'start' component will determine when 
                          ## harm monitoring beings.
                          harmMonitorRange, 
                  
                          ## Total Type I error for potential harm monitoring (per vacc. arm)
                          harmMonitorAlpha=0.05,
                  
                          ## 'harmMonitorContol' is a list containing other parameters controlling 
                          ##    the harm monitoring. Currently the only one used in 'maxCnt'.  This
                          ##    *MUST* be specified and represesents the set of infection counts for
                          ##    which harm bounds will be constructed.  Harm monitoring ends when
                          ##    'maxCnt' has been exceeded, whether or not you've reached the
                          ##    criteria for initiaton of non-efficacy monitoring, so you should: 
                          ##    (a) choose a nonEff start method that guarantees starting
                          ##        at/before 'maxCnt', or
                          ##    (b) choose maxCnt to be a value larger than you're non-efficacy
                          ##        starting count will ever be (maybe use 2 or 3 times the 
                          ##        upper range for harm monitoring (harmMonitorRange[2]).
                          ##  If maxCnt is left unspecified, the default value of 
                          ##    3 x harmMonitorRange[2] will be used for it.
                          #harmMonitorControl=list( maxCnt=NULL),
                  
                          ## if you have determined the constant alpha value to do the binomial
                          ## 'potential-harm-monitoring' tests at, they you can pass it through
                          ## this argument, in which case argument 'harmMonitorAlpha' need not
                          ## be specified
                          alphaPerTest=NULL,
                  
                          ## character vector of methods to use to determine when to start non-efficacy
                          ## monitoring.  Each method requires different input information, and that
                          ## information should be input via the 'nonEffStartParams' argument - which
                          ## should be a list that the specifies the inputs for the method.  Each
                          ## method's input requirements should be specified in the documentation.
                          ##
                          ## Methods:
                          ##  FKG - the starting method suggested in Freidlin, Korn and Gray's 2010
                          ##        'Clinical Trials' paper:  "A general inefficacy interim 
                          ##        monitoring rule for randomized clinical trials"
                          ##        It boils down to start monitoring at the earlier infection count
                          ##        such that an estimated effect <= 0 would cause the trial to
                          ##        stop.  In our code that translates to a 95% CI (based on the
                          ##        asymptotic variance of the log-rank statistic) around an 
                          ##        estimated VE of 0% would exclude the VE specified by parameter
                          ##        'upperVEnonEff'.  The expectation is that this will be harmonized
                          ##        with the non-efficacy monitoring, and so will use the arguments
                          ##        'upperVEnonEff' 'alphaNoneff' from argument list.  However, should
                          ##        you want to use different values, you may do say by passing them
                          ##        via the 'nonEffStartParams' argument list.  The same naming should
                          ##        be used within the list.
                          ##      Parameters: upperVEnonEff alphaNoneff  
                          ##  
                          ##  fixed - Starts at the infection count specified by paramter 'N1'. 
                          ##          If 'N1' is set to 75, then non-efficacy always begins at the
                          ##          75th infection 
                          ##      Parameter: N1
                          ##
                          ##  ? - method specifies start time by a combination of variables:
                          ##      Parameters:
                          ##         minCnt - gives the minimum (combined) infection count at which
                          ##               monitoring can being
                          ##         maxCnt - gives the maximum (combined) infection count at which
                          ##               monitoring can being (it will begin at this count 
                          ##               whether or not the other criteria are satisfied).
                          ##               [OPTIONAL]
                          ##         lagTimes, lagMinCnts - these arguments are paired, each can
                          ##               be a vector of length > 1. 'lagMinCnts' is a vector
                          ##               if infection counts that are required to occur 'later
                          ##               in the trial' (larger value of follow-up time) than
                          ##               the associated time given in vector 'lagTimes'.
                          ##               E.g if lagTimes=c(26,39) and lagMinCnts=c(25,5) this
                          ##               says we cannot start until we've reached 25 infections
                          ##               that have occurred *after* (not including) week 26, and
                          ##               also 5 that have occurred after week 39.
                          ##       NOTE: Parameters: minCnt, lagTimes and lagMinCnts are all required,
                          ##               (minCnt can be set to 0 if desired); maxCnt is optional
                          ##
                          ##
                          ##  custom - [ NOT YET IMPLEMENTED ] You provide the function and the parameters
                          ##      Parameters: 'Func' (your function) + whatever parameters your function needs
                          ##
                          ##  old - the old 'seqDesign' starting method.  Severely deprecated.
                          ##      Parameters: minCnt, minPct, week1, week2
                          ##
                          nonEffStartMethod=c("FKG", "fixed", "?", "old"),
                  
                          ## A *list* (not a vector) of parameters needed by the method specified
                          ## in argument 'nonEffStartMethod'.  Some methods have defaults in place,
                          ## for those you do not need to use 'nonEffStartParams' unless you want 
                          ## values other than the defaults.
                          nonEffStartParams=NULL,
                  
                          #minCnt,
                          #maxCnt, ## the maximum start of non-efficacy monitoring
                          #minPct,
                          #week1,
                          #minCnt2,
                          #week2,
                  
                          nonEffIntervalUnit=c("counts","time"),
                          nonEffInterval,
                  
                          ## lowerVEnoneff is not required.  Specify only if you want this
                          ## condition as part of your monitoring.
                          lowerVEnoneff=NULL,
                          upperVEnoneff,
                          highVE, 
                          stage1VE,
                          lowerVEuncPower=NULL,  
                  
                          alphaNoneff,
                          alphaHigh,
                          alphaStage1,
                          alphaUncPower=NULL,
                          
                          estimand=c("combined", "cox", "cuminc"),
                  
                          ## 'laggedMonitoring' replaces argument 'post6moMonitor' and 
                          ## 'lagTime' replaces 'VEcutoffWeek'
                          laggedMonitoring=FALSE,
                          lagTime = NULL,
                  
                          saveFile= NULL,
                          saveDir = NULL,
                          verbose = TRUE ) {
  ## selected Arguments
  ##
  ##   lowerVEnoneff: 
  ##       lower confidence bound to be below 'lowerVEnoneff' to meet criterion 1
  ##       for declaring non-efficacy
  ##
  ##   upperVEnoneff:
  ##       upper confidence bound to be below 'upperVEnoneff' to meet criterion 2
  ##       for declaring non-efficacy
  ##
  ##   highVE:
  ##       lower confidence bound to be above 'highVE' for declaring high efficacy
  ##
  ##   stage1VE:
  ##       lower confidence bound to be above 'stage1VE' for advancing into Stage 2
  ##
  ##   lowerVEuncPower: 
  ##       vector of VE value for which the user wishes to determine unconditional
  ##       power to reject  H0: VE <= (lowerVEuncPower * 100)% 
  ##       at 2-sided alpha level 'alphaUncPower'
  ##
  ##   alphaNoneff:
  ##       One minus confidence level of 2-sided CI for non-efficacy monitoring
  ##
  ##   alphaHigh:
  ##       One minus confidence level of 2-sided CI for high efficacy monitoring
  ##
  ##   alphaStage1:
  ##       One minus confidence level of 2-sided CI for testing whether an arm 
  ##       advances into Stage 2
  ##
  ##   alphaUncPower:  
  ##       One minus confidence level of 2-sided CI for unconditional power to
  ##       reject H0: VE <= (lowerVEuncPower * 100)% 

  # for reverse compatibility
  post6moMonitor <- laggedMonitoring
  VEcutoffWeek <- lagTime

  estimand <- match.arg(estimand)
  nonEffIntervalUnit <- match.arg(nonEffIntervalUnit)
  nonEffStartMethod  <- match.arg(nonEffStartMethod)
  
  if ( is.list(dataFile) ) {
    trialObj <- dataFile
  } else {
    if ( !is.null(saveDir) ){
      ## load in RData object (a list named 'trialObj' )
      load(file.path(saveDir, dataFile))
    } else {
      load(dataFile)
    }
  }

  ## check contents of 'trialData'
  if ( !all( c("entry","exit","event","trt") %in% names(trialObj$trialData[[1]]) ) )
    stop("Trial data must contain columns: entry, exit, event, trt\n")
  
  nTrtArms <- as.integer( trialObj$nArms - 1 )
  nTrials <- length(trialObj[["trialData"]])

  ## set altVE based on lowerVEuncPower, etc. if those were specified instead
  alphaAltVE <- altVE <- NULL
  if ( is.null(altVE) && !is.null(lowerVEuncPower)) {
      altVE <- lowerVEuncPower
  }    
  if (is.null(alphaAltVE) && !is.null(altVE)) {
      alphaAltVE <- ifelse(!is.null(alphaUncPower), alphaUncPower, alphaStage1)
  }

  ## ----------------------------------------------------------------------
  ##   The following code is needed to determine infection count at which to 
  ##   begin the nonEff monitoring.  Also needed for harm-monitoring

  ## derive the prob. of assignment to each vaccine arm (relative to placebo
  ## only, not relative to other vaccine arms)
  null.p <- numeric(length=nTrtArms)
  for (i in 1:nTrtArms) {
      arms.i <- c(1, i+1)
      assn.probs <- trialObj$trtAssgnProbs[ arms.i ]
      null.p[i] <- assn.probs[2] / sum(assn.probs)
  }

  ## if we have multiple vaccine groups and they are of different sizes,
  ## then we must stop.  The  code is not set up to handle that yet.
  if (nTrtArms > 1  && ( max(null.p)-min(null.p) > sqrt(.Machine$double.eps) ) )
  {
     stop("The code is not currently set up to handle harm monitoring when\n",
          "there are multiple treatment arms and their sizes differ. \n\n")
  } else {
    ## all values of null.p are equal, just choose the first
    null.p <- null.p[1]
  }


  ## This will also determine when nonEfficacy monitoring starts
  ## (and therefore how long harm monitoring will last).
  ## Two methods are found here and two later on.  These two do not
  ## require any specific information about the trial to determine,
  ## whereas the later methods do.

  ## Method implemented here:
  ##    Start nonEff monitoring at the smallest infection count for which we 
  ##    would be able to reject the null hypothesis that VE = upperVEnoneff
  ##    if our point estimate of VE = 0.  This computation requires us to 
  ##    assume that our test statistic is asymptotically equivalent to the
  ##    log-rank statistic, which (under the null of no diff - not the setting
  ##    we're in, right?) has a distribution that is 
  ##      Normal( 0, sqrt(1/(p*(1-p)*D))),  where
  ##    D = number of endpoints between the placebo and trt arm being monitored,
  ##    p = prob. of assignment to the vaccine arm vs. the placebo arm (all 
  ##        other arms ignored for computing 'p')
  ##
  if ( nonEffStartMethod == "FKG") {
 
     if (!is.null(nonEffStartParams)) {  
        ## values specified directly by user (maybe - we don't
        ## know what's actually in 'nonEffStartMethod' yet
        upperVE  <- nonEffStartParams$upperVEnoneff
        alphaNE  <- nonEffStartParams$alphaNoneff

        if ( is.null(upperVE) && is.null(alphaNE) )
            warning("The argument 'nonEffStartParams' has been specified but",
                    "it does not contain\n", "parameters used by method FKG\n\n",
                    immediate.=TRUE)
     } else {
        upperVE <- NULL
        alphaNE <- NULL
     }

     ## Fill in NULL values for either parameter using the monitorTrial 
     ## arguments upperVEnoneff and alphaNoneff
     if ( is.null(upperVE) )
         upperVE <- upperVEnoneff

     if ( is.null(alphaNE) )
         alphaNE <- alphaNoneff

     if ( is.null(upperVE) || is.null(alphaNE) )
         stop( "One or both of the arguments needed by nonEffStartMethod 'FKG'",
            "are missing\n.", "They were not specified via argument ",
            "'nonEffStartParams' and their values also could not be obtained\n",
            "from the argument 'upperVEnoneff' and 'alphaNoneff'\n",
            "Exiting...\n" )
 
     N1 <- ceiling( qnorm(1 - alphaNE/2)^2 /
                     ( null.p*(1-null.p)*log(1-upperVE)^2 ) )

     ## set upperbound of 'harmMonitorRange' to be 'N1'
     harmMonitorRange <- c(harmMonitorRange[1], N1)    

  } else if ( nonEffStartMethod == "fixed" ) {

     if (!is.null( nonEffStartParams$N1 ) ) {
        N1 <- nonEffStartParams$N1

        ## set upperbound of 'harmMonitorRange' to be 'N1'
        harmMonitorRange <- c(harmMonitorRange[1], N1)    
     } else {
        stop("Argument 'nonEffStartMethod' was specified as 'fixed'. ",
             "This method requires\n", "you to provide an argument named",
             "'N1' via the list argument 'nonEffStartParams'\n.",
             "e.g. nonEffStartParams <- list(N1=60).  Please fix.\n" )
     }
  }


  ## -------------------------- get harm bounds ---------------------------  

  # calculate stopping boundaries for harm
  ## NOTE: harmMonitorRange dictates that range over which the type-I error is
  ## spent, it does not specify the time at which harm monitoring stops.  That is
  ## dictated by when non-eff monitoring begins
  if ( is.null(alphaPerTest) ){ 
    ## choose the value of alphaPerTest to spend the type I error over harmMonitorRange
    alphaPerTest <- getAlphaPerTest(harmMonitorRange, null.p,
                                    totalAlpha = harmMonitorAlpha)
  }

  ## If nonEffStartMethod "?" was chose then see if maxCnt was used
  if (nonEffStartMethod == "?" && !is.null(nonEffStartParams$maxCnt) ) {
      if ( is.infinite(nonEffStartParams$maxCnt) ) 
          nonEffStartParams$maxCnt <- NULL
      else
          maxCnt <- nonEffStartParams$maxCnt
  }

  ## For generating the bounds, we need both 'N' and the upper limit of 
  ## 'harmMonitorRange' to be an upper limit on the number of tests to do
  ## (i.e. an upper limit on when the first non-efficacy will be done)
  if (exists("N1")) {
    maxCnt = N1

  } else if ( !exists("maxCnt") ) {
    maxCnt = 5 * harmMonitorRange[2]
  }

  ## get harm bounds 
  harmBounds <- getHarmBound(
                    N = maxCnt,
                    per.test = alphaPerTest, 
                    harmBoundRange = c(harmMonitorRange[1], maxCnt),
                    null.p = null.p, 
                    dataDir = saveDir, 
                    verbose = verbose)

  ## -------------------------- end get harm bounds ---------------------------  


  ## creates a list of length 'nTrials' each element of which is a list of 
  ## length 'nTrtArms'
  out <- rep( list(vector("list",nTrtArms)), nTrials )


  ## censor all trials to stage 1 (fastest to do all at once).  Still need to
  ## keep access to uncensored data, for highEff monitoring
  if ( !missing(stage1) ) {
      stg1dat <- censorTrial(trialObj$trialData, times=stage1,
                             timeScale="follow-up")
  }

  ## Start looping over trials
  for (i in 1:nTrials ) {

    ## extract data for the i-th trial
    if ( !missing(stage1) ) {
       datI <- stg1dat[[i]]
       datIall <- trialObj[["trialData"]][[ i ]]
    } else {
       datI <- trialObj[["trialData"]][[ i ]]
    }
    
    ## create separate set of only infections ('events'), then
    ## order the events by trial time at which they are observed
    eventDF <- datI[datI$event == 1, ]
    eventDF <- eventDF[order(eventDF$exit), ]

    #cat("Trial :", i, "  \n")

    ## Now start comparisons - each active arm vs. placebo arm, one at a time
    for (j in 1:nTrtArms) {
      
      ## subset *all data* down to the two arms being compared
      datI.j <- datI[datI$trt %in% c(0,j), ]
      
      ## convert 'trt' to indicator variable before passing it to
      ## 'applyStopRules' (i.e. convert the non-zero values to 1)
      datI.j$trt <- as.integer(datI.j$trt > 0 )

      ## same for the 'all' data.frame
      if ( exists("datIall") ){
         datIall.j <- datIall[datIall$trt %in% c(0,j), ] 
         datIall.j$trt <- as.integer(datIall.j$trt > 0 )
      }
      


      ## 1. do harm monitoring for j-th treatment arm vs. placebo

      ## subset events in relevant arms: j-th active trt and placebo(trt=0)
      E.j <- eventDF[eventDF$trt %in% c(0,j), ]
      nInfec <- nrow(E.j) # counts infections through 'stage1'
      E.j$nInf <- 1:nInfec

      if (nonEffStartMethod == "?") {
          if (!is.null(nonEffStartParams)) {
              argNames <- c("minCnt","maxCnt","lagTimes","lagMinCnts")
              argList <- nonEffStartParams[ 
                             argNames[ argNames %in% names(nonEffStartParams) ] ]

              if (!all( c("minCnt","lagTimes","lagMinCnts") %in% names(argList) ))
                  stop("The argument 'nonEffStartParams' has been specified but",
                       "it does not contain\n", "all the parameters needed by ",
                       "nonEffStartMethod '?'\n") 
          } else {
              stop("The use of nonEffStartMethod ", nonEffStartMethod, " requires",
                   " that you specify parameters: \n",
                   " 'minCnt', 'lagTimes' and 'lagMinCnts' \n",
                   "They must be passed to 'monitorTrial' through the argument ",
                   "nonEffStartParams - which must be a list.\n",
                   "Please fix and then rerun\n")
          }
          ## gets infection count at which non-efficacy monitoring begins
          N1 <- getFirstNonEffCnt(datI.j, minCnt=argList$minCnt,
                    lagTimes=argList$lagTimes, lagMinCnts=argList$lagMinCnts )

      } else if (nonEffStartMethod == "old") {

          if (!is.null(nonEffStartParams)) {
              ## values specified directly by user - maybe. 
              ## We don't know what's actually in 'nonEffStartMethod' yet
              minPct  <- nonEffStartParams$minPct
              minCnt  <- nonEffStartParams$minCnt
              week1   <- nonEffStartParams$week1
              minCnt2 <- nonEffStartParams$minCnt2
              week2   <- nonEffStartParams$week2

              if ( is.null(minPct)  || is.null(minCnt) || is.null(week1) ||
                   is.null(minCnt2) || is.null(week2)  ) 
                  stop("The argument 'nonEffStartParams' has been specified but",
                      "it does not contain\n", "all the parameters needed by ",
                      "nonEffStartMethod 'old'\n")
          }

          ## Old method of determining "N1" - depracated for mulitple reasons:
          ## it cannot easily be determined when/if the minPct criteria will be
          ## satisfied, introducing an inability to effectively anticipate and
          ## plan for analyses.  Other methods allow for use of various forms 
          ## of minimum infection cutoffs - those should be used instead
          N1 <- getInfecCntFirstNonEff( E.j, 
                     minPercent = minPct,
                     minCount = minCnt, 
                     week1 = week1, 
                     nInfecAfterwk = minCnt2,
                     week2 = week2)
      }

      ## Ensure we have harm bounds up to value 'N1', if not, regenerate them 
      ## so we have enough.
      if (N1 > max(harmBounds$N)) {
          maxCnt <- 2 * N1

          ## get harm bounds 
          harmBounds <- getHarmBound(
                            N = maxCnt,
                            per.test = alphaPerTest,
                            harmBoundRange = c(harmMonitorRange[1], maxCnt),
                            null.p = null.p,
                            dataDir = saveDir,
                            verbose = verbose)
      }

      ## subset to restrict harm monitoring to the first 'N1' infections
      harmBounds.j <- harmBounds[harmBounds$N <= N1, c("N","V","P")]
      
      ## doEvaluate harm monitoring bounds.  Just based on infected ppts, so
      ## only uses 'E.j'
      harmRes  <- do_harm_monitoring(E.j, bounds = harmBounds.j ) 
      
      ## if there is indication of "harm" then need to create store output
      ## object then move to next trial 
      if ( harmRes$isHarm ) {
          fst <- finalStage1Test(
                     datI.j,
                     analysisType = "stopTime",
                     lowerVE = altVE,
                     alphaLevel = alphaAltVE,
                     estimand = estimand,
                     time = harmRes$stopTime,
                     randFraction = null.p )

          fest <- fst$stage1summ 
          VEcir <- !is.null(fest$VE_CIR)
          VE_CI <- if (VEcir) {
                       c(fest$VE_CIR_loCI, fest$VE_CIR_upCI)
                   } else {
                       c(fest$VE_Cox_loCI, fest$VE_Cox_upCI)
                   }

          out[[i]][[j]] <- 
              list( 
                  boundHit = "Harm",
                  stopTime = harmRes$stopTime, 
                  stopInfectCnt = harmRes$stopInfectCnt,
                  stopInfectSplit = harmRes$stopInfectSplit,
                  stopVE = ifelse(VEcir, fest$VE_CIR, fest$VE_Cox),
                  stopVE_CI= VE_CI,
                  VE_estimand = estimand,
                  stage1complete = FALSE,
                  stage2complete = FALSE,
                  altDetected = fst$altDetected  )

        next
      } 


      ## If harm bounds not hit, move to nonEff monitoring

      ## 2. begin non-eff monitoring prep
      ## --------------------------------


      if (N1 > nInfec) {

          ## deal with annoying case where the trial didn't have enough infec.s
          ## to reach the start point of the non-eff. interim analyses
           
          ## For now we do nothing, except skip all the interim analysis code
          ## and go directly to the finalStage1 test
          nonEffTimes <- NA

          ## create an empty object futRes needed later to make code work
          futRes <- list()

          ## Set 'N1' to NA ??  Need to think about it
          ## N1 <- NA

      } else {  # This 'else' goes on FOREVER!
    
          ## determine 'nonEffTimes' and then do non-eff monitoring
          ## Determine the times at which nonEff monitoring will occur
          if (nonEffIntervalUnit == "counts") {

            ## makes sequence of counts at which analyses will be done
            if (length(nonEffInterval)==1){
              # then a constant increment in the endpoint count is assumed
              nonEffCnts <- seq(from = N1, to = nInfec, by = nonEffInterval)  
            } else {
              # then 'nonEffInterval' is a vector of increments specifying endpoint counts for subsequent interim looks
              nonEffCnts <- c(N1, N1 + nonEffInterval)
            }

            ## Convert counts into times.  
            nonEffTimes <- getInfectionTimes(datI.j, cnts=nonEffCnts )
          } else {
              ## User specified a time-sequence for monitoring. Figure out when 
              ## the first time will be, and the time of the last infection too,
              ## then create the sequence
              firstLastnonEffTimes <- 
                  getInfectionTimes(datI.j, cnts=c(N1, nInfec))

              nonEffTimes <- 
                  seq(from = firstLastnonEffTimes[1],
                        to = firstLastnonEffTimes[2],
                        by = nonEffInterval)
          }

          ## run non-eff 
          futRes <- 
              applyStopRules (
                  datI.j,
                  testTimes = nonEffTimes,
                  boundType = "nonEff",
                  boundLabel = "NonEffInterim", 
                  lowerVE = lowerVEnoneff,
                  upperVE = upperVEnoneff,
                  alphaLevel = alphaNoneff,
                  laggedMonitoring = laggedMonitoring,
                  lagTime = lagTime,
                  estimand=estimand,
                  randFraction = null.p )

          ## if noneff hit, then store results and go to next arm
          ## (note: output will be a list because 'futRes' is a list
          if ( futRes$boundWasHit ) {

             fst <- finalStage1Test(
                        datI.j,
                        analysisType = "stopTime",
                        lowerVE = altVE,
                        alphaLevel = alphaAltVE,
                        estimand = estimand,
                        time = futRes$stopTime,
                        randFraction = null.p )

              fest <- fst$stage1summ 
              VEcir <- !is.null(fest$VE_CIR)
              VE_CI <- if (VEcir) {
                           c(fest$VE_CIR_loCI, fest$VE_CIR_upCI)
                       } else {
                           c(fest$VE_Cox_loCI, fest$VE_Cox_upCI)
                       }

              out[[i]][[j]] <- 
                  c( futRes,
                     list( stopVE = ifelse(VEcir, fest$VE_CIR, fest$VE_Cox),
                           stopVE_CI= VE_CI,
                           VE_estimand = estimand,
                           firstNonEffCnt = N1,
                           stage1complete = FALSE,
                           stage2complete = FALSE,
                           altDetected = fst$altDetected ) )
              next
          } 


          ## 3. if nonEff NOT hit, then move on to high-eff monitoring
          ##    We first do *only* at the stage1 times (i.e. same times as for
          ##    nonEff monitoring ), but using all data  
          highEffRes <- applyStopRules(
                            datIall.j,
                            testTimes = nonEffTimes,
                            boundType = "highEff",
                            boundLabel = "HighEff", ## assumed by later function
                            lowerVE = highVE,
                            alphaLevel = alphaHigh,
                            laggedMonitoring = FALSE,
                            estimand=estimand,
                            randFraction = null.p )

          if ( highEffRes$boundWasHit ) {
              fst <- finalStage1Test(
                         datI.j, 
                         analysisType = "stopTime", 
                         lowerVE = altVE, 
                         alphaLevel = alphaAltVE, 
                         estimand = estimand, 
                         time = highEffRes$stopTime,
                         randFraction = null.p )

              fest <- fst$stage1summ 
              VEcir <- !is.null(fest$VE_CIR)
              VE_CI <- if (VEcir) {
                           c(fest$VE_CIR_loCI, fest$VE_CIR_upCI)
                       } else {
                           c(fest$VE_Cox_loCI, fest$VE_Cox_upCI)
                   }

              ## store output, then determine value for 'altDetected'
              ## (note: output will be a list because 'highEffRes' is a list
              out[[i]][[j]] <- 
                  c( highEffRes,
                     list( summNonEff = futRes$summObj,
                           stopVE = ifelse(VEcir, fest$VE_CIR, fest$VE_Cox),
                           stopVE_CI= VE_CI,
                           VE_estimand = estimand,
                           firstNonEffCnt = N1,
                           stage1complete = FALSE, 
                           stage2complete = FALSE,
                           altDetected = fst$altDetected ) )

              next
          }
      }  ## end of (N1 > nInfec) if/else clause

##-----------------------------------------------------------------------

      ## Reached the end of stage 1 w/o hitting any boundaries (yee-ha!)
      ## celebrate by doing 'end-of-stage-1' dance 
      fst <- finalStage1Test(
                 datI.j, 
                 analysisType = "final", 
                 stage1VE = stage1VE,
                 lowerVE = altVE,
                 alphaLevel = alphaStage1, 
                 estimand = estimand, 
                 boundLabels=c("Eff", "NonEffFinal"),
                 randFraction = null.p )

      fest <- fst$stage1summ 
      VEcir <- !is.null(fest$VE_CIR)
      VE_CI <- if (VEcir) {
                   c(fest$VE_CIR_loCI, fest$VE_CIR_upCI)
               } else {
                   c(fest$VE_Cox_loCI, fest$VE_Cox_upCI)
               }

      ## not storing into out[[i]][[j]] yet, we'll add on if not final-non-eff
      finalStg1List <- 
          c( fst,
             list( stopVE = ifelse(VEcir, fest$VE_CIR, fest$VE_Cox),
                   stopVE_CI= VE_CI,
                   VE_estimand = estimand,
                   firstNonEffCnt = N1,
                   summNonEff  = futRes$summObj,
                   stage1complete = TRUE,
                   stage2complete = FALSE ) )

      ## if we hit the bound at stage 1 analysis, then we're done with this cohort
      if (fst$boundHit == "NonEffFinal") {
          out[[i]][[j]] <- finalStg1List
          next
      }

      ## 5.  If no stage1 highEff, and we "passed" the end-of-stg1 analysis, then do
      ##    the stage2 highEff analysis

      ## (a) timing of it is  half way between end of stage 1 and end of trial 
      ##     compute as end-of-stg1 (max(datI.j$exit)) plus half the remaining time
      ##     (stage2- stage1)/2 until the last person from stage1 completes stg2
      stg2highEffTime <- ( max(datI.j$exit) + (stage2 - stage1)/2 )

      stg2highEff <- 
          applyStopRules(
                 datIall.j,
                 testTimes = stg2highEffTime,
                 boundType = "highEff",
                 boundLabel = "HighEff", ## assumed by later function
                 lowerVE = highVE,
                 alphaLevel = alphaHigh,
                 laggedMonitoring = FALSE,
                 estimand=estimand,
                 randFraction = null.p )


      ## Note: 'altDetected' result is included in 'finalStg1List', which is why
      ##  you will not see it added to any of the following output objects
      if ( stg2highEff$boundWasHit ) {
          
          ## alter some values in present in 'finalStg1List'
          finalStg1List$boundHit  <- stg2highEff$boundHit
          finalStg1List$boundType <- stg2highEff$boundType
          finalStg1List$stopTime  <- stg2highEff$stopTime

          ## remove components from stg2Eff already present in 'finalStg1List'
          dropComps <- c("boundWasHit", "boundHit", "boundType", "stopTime")
          stg2HE    <- stg2highEff[ !names(stg2highEff) %in% dropComps ]

          ## Note: 'altDetected' result is included in 'finalStg1List'
          out[[i]][[j]] <- c( finalStg1List, 
                              list( summNonEff = futRes$summObj ),
                              stg2HE ) 
      } else {
          ## Stage 2 completed without hitting high-efficacy
          finalStg1List$stage2complete <- TRUE
          finalStg1List$stopTime       <- max( datIall.j$exit)
          out[[i]][[j]] <- c( finalStg1List,
                              list( 
                                summNonEff  = futRes$summObj,
                                stg2highEff = stg2highEff$summObj ) )
      }
    }
  }
##-----------------------------------------------------------------------
#   bound            stage1complete  altDetected(stage1 tests - 2 of them)
#   "NonEffInterim"  FALSE           FALSE 
#   "HighEff"        FALSE           TRUE/FALSE [could be either]
#   "NonEffFinal"    TRUE            FALSE 
#   "Eff"            TRUE            TRUE/FALSE [could be either]
#   "HighEff"        TRUE            (Same as for "Eff") 
##-----------------------------------------------------------------------
  
  ## I left this in here, but not very useful - will excise later probably
  if (verbose){
    for (i in 1:nTrtArms) {
      cat("Probabilities of reaching each possible conclusion:\n")
      print( round( table(sapply( out, function(x) x[[i]]$boundHit ),
                          useNA="ifany")/nTrials, 4))
      
      if (!is.null(altVE)){
        altDetectedMatrix <- 
            do.call("rbind", 
                    lapply(out, function(x) x[[i]]$altDetected ))

        # this occurs if each trial stops for potential harm, and thus the value of 
        # 'altDetected' is NULL for each trial
        if ( is.null( altDetectedMatrix ) ) {
            designPower <- 0
        } else {
            designPower <- colSums(altDetectedMatrix, na.rm=TRUE)/nTrials
        }
        cat("\nUnconditional power to reject the specified null hypotheses =",
             format(designPower, digits=3, nsmall=3), "\n\n")
      }
    }
  }
  
  ## save monitoring output
  if ( !is.null(saveDir) ) {
    if ( is.null(saveFile) ) {
        if ( is.list(dataFile) )
          warning(
              "The output of 'monitorTrial' will not be saved to a file.\n",
              "You have not specified the argument 'saveFile' and a default\n",
              "filename cannot be constructed when argument 'dataFile' is ",
              "a list.\n\n", immediate.=TRUE)
        saveFile <- paste0("monitorTrial", 
                           substr(dataFile, 9, nchar(dataFile)-6), 
                           "_", estimand, ".RData")
    }
    save(out, file = file.path(saveDir, saveFile) )
    if (verbose) { 
        cat("Output saved in:\n", file.path(saveDir, saveFile), "\n\n") 
    }
  } 

  ## it should not be an "either/or" decision whether you save or output the results
  return( invisible( out ) )
}

######################################## end of 'monitorTrial'#########################################




### This version is Sept 10, 2015
###   I'm beginning major modification of the function to eliminate stuff that
###   is unnecessary and makes no sense.
###     - only one set up "upper","lower" parameters to be used for both 
###       high eff and noneff.  Only one can be done at a time anyway.
###     - the upper/lower arg.s will be in terms of VE and not HR, since all of
###       thinking is done on the VE scale (and HR only applies to cox model)
###     - removing all 'final analysis' related things, which will include
###       the 'UncPower' params, Stage1 params, etc.
###     - only one 'alpha' parameter
###     - a more extensive set of output information
###     - < etc >

applyStopRules <- 
    function(d, infectionTotals=NULL, testTimes=NULL, 
             boundType = c("nonEff","highEff"), boundLabel=boundType,
             lowerVE=NULL,  upperVE=NULL, alphaLevel=0.05,
             laggedMonitoring=FALSE,  lagTime=NULL, laggedOnly=FALSE,
             estimand = c("cox","cuminc","combined"),
             randFraction)
{
  
  ## This function apply the stopping rules for nonefficacy and high efficacy monitoring.
  ##
  ## Argument 'd' should be a data.frame containing (at least) columns:
  ##   'entry' - the entry time (in trial time)
  ##   'exit'  - time of event, trial completion or dropout (trial time)
  ##   'event' - 0/1 indicator of event of interest
  ##   'trt'   - 0/1 indicator of vaccine receipt
  ## 
  ## Other arguments:
  ## ---------------
  ## infectionTotals - a vector specifying the total number of infections
  ##                     at which analyses should take place.  
  ## 
  ## testTimes   - a vector specifying the calendar times at which noneff analyses
  ##                 should take place.  Either this arg. or 'infectionTotals' must
  ##                 be specified.  This arg is preferred as it will need to be
  ##                 deteremined within the function if not given by the user.
  ## 
  ## boundType   - character string indicating the type of bound being evaluated
  ##               must be one of "nonEff" or "highEff".
  ##
  ## boundLabel  - a label to use for naming the bound.  Can be whatever you like.
  ##               Defaults to the value of argument 'boundType'
  ##
  ## lowerVE     - if boundType=='nonEff' then this specifies the VE value that the
  ##               lower CI of the estimated VE(s) must lie below for non-efficacy
  ##               stopping to occur.   If boundType='highEff' then this specifies
  ##               the VE value that the lower CI of the estimated VE(s) must lie 
  ##               ABOVE for high-efficacy stopping to occur
  ## 
  ## upperVE     - applicable only when boundType=="nonEff".  Specifies the VE that
  ##               the upper CI of the estimated VE(s) must lie below in order to
  ##               stop for nonEfficacy.
  ## 
  ## alphaLevel  - two-sided alpha level to use in constructing 95% confidence 
  ##               interval for the VE(s) estimated by the estimand(s) specified
  ##               by the argument 'estimand'
  ##
  ## laggedMonitoring - TRUE/FALSE (default is FALSE)
  ##               Should monitoring be done on a subset of data that excludes
  ##               the first 'lagTime' weeks of follow-up for each participant.
  ##
  ## lagTime     - The amount of "lag-time" before lagged-monitoring begins.
  ##               Must be specified if laggedMonitoring is TRUE
  ##
  ## laggedOnly - TRUE/FALSE, should *only* lagged monitoring be done            
  ##               (not currently implemented)
  ##
  ## estimand    - a character string specifying the estimand(s) to be used in
  ##                 monitoring (can be one of "combined", "cox", and "cuminc")
  ##
  ## randFraction - the fraction of randomizations going to a particular 
  ##                treatment arm, when only that treatment arm and the placebo 
  ##                arm are considered.
  ## ---------------------------------------------------------------------------

  ## check contents of 'd'
  if ( !all( c("entry","exit","event","trt") %in% names(d) ) )
    stop("DataFrame 'd' must contain columns: entry, exit, event, and trt\n")

  boundType <- match.arg(boundType)
  estimand  <- match.arg(estimand)
  
  if (boundType == "nonEff"){ 
    if ( any( is.null(upperVE), is.null(alphaLevel) ) ) {
        stop("The arguments 'upperVE' and 'alphaLevel' must be specified for",
             " non-efficacy monitoring.\n") 
     }
  } else {
    # high-efficacy monitoring
    if ( any( is.null(lowerVE), is.null(alphaLevel) ) ) {
        stop("The arguments 'lowerVE' and 'alphaLevel' must be specified for ",
             "high efficacy monitoring.\n")
    }
  }

  ## Make sure we have only two trt groups and that they're coded as 0 and 1
  uniq.trt <- sort( unique(d$trt) )
  if ( length(uniq.trt)>2 ) {
    warning("The data set given to 'applyStopRules' contains", length(uniq.trt),
            "treatments - it should only contain 2.\n\n", immediate.=TRUE) 
  } else if ( any( uniq.trt != c(0,1) ) ) {
    warning("The data set given to 'applyStopRules' contains values of 'trt'",
            "other than 0 and 1, this is probably an error.\n",
            "The values are: ", uniq.trt[1], " and ", uniq.trt[2], "\n\n",
             immediate.=TRUE)
  }

  ## If 'lowerVE' is NULL we set it to a large value so that we can use a common
  ## set of code (no special cases) and not have the bound have any effect - 
  ## since the associated criteria will always be met. Note: max possible VE is 1

  if ( is.null(lowerVE) ) 
    lowerVE <- 17

  ## force evaluation of some input arguments to avoid scoping issues when passing
  ## them down into subfunctions.
  force( alphaLevel )
  force( upperVE )
  force( randFraction )


  ### -------- Stop checking, start working ---------



  ## If 'infectionTotals' was given, then need to get associated times
  if ( is.null(testTimes) ) {
    if ( !is.null(infectionTotals) ) {
        testTimes <- getInfectionTimes(d, cnts=infectionTotals)  
    } else {
      stop("One of the arguments ('infectionTotals', 'testTimes') must be specified\n") 
    }
  }

  ## the length of testTimes indicates how many tests we'll be doing
  nTests <- length(testTimes)

  if (estimand == "combined" && boundType=="highEff") {
      estimand <- "cuminc"
  }

  ## set indicators for estimation types to use
  cox  <- (estimand %in% c("cox",   "combined"))
  cir  <- (estimand %in% c("cuminc","combined"))


  ## create an indicator vector of what is being evaluated in the monitoring
  ## (i.e. which estimates and which data - non-lagged/lagged )
  indEst <- c(cir = cir,  
              cox = cox, 
              lagcir = cir && laggedMonitoring,
              lagcox = cox && laggedMonitoring )

  ## define function for testing stopping criteria based on boundType
  ## note: the as.vector() is included to strip off any names attached
  evalStopCrit <- if ( boundType=="nonEff" ) {
                    function(VEci, lowerVE, upperVE) 
                      return( as.vector( (VEci[,1] < lowerVE) & (VEci[,2] < upperVE)) )
                  } else {
                    ## boundType == "highEff"
                    function(VEci, lowerVE, ...) 
                      return( as.vector( VEci[,1] > lowerVE ) )
                  }

  ## set to prevent R from changing data.frame character columns into factors
  options( stringsAsFactors = FALSE)

  ## model formula to be used by coxph() and/or survfit()
  survFormula <- Surv(exit-entry, event) ~ trt
  

  ## ***** Create list(s) containing all datasets needed for nonefficacy *******
  
  ## (1) Create a list of censored datasets - one per element of 'testTimes'
  censDatList <- censorTrial(d, times=testTimes, timeScale="calendar")

  ## (2) If laggedMonitoring == TRUE, create a 2nd list of censored datasets,
  ##     this time left censored to 'lagTime' (based on follow-up time).  Note
  ##     that this set is based off the right censored data in censDatList.
  if (laggedMonitoring) {
      censDatList_lag <- censorTrial(censDatList, times=lagTime, 
                                          timeScale="follow-up", type="left")
  }


  ## ***** Compute estimators needed for requested monitoring  *****

  ## cumulative incidence ratio, if requested
  if (cir) {
      cumIncOut <- cumIncRatio( censDatList, times="last", 
                                alphaLevel = alphaLevel,
                                randFraction = randFraction)

      CIRobj <- transform( cumIncOut, 
                    infectTotal = nEvents,
                    infectSplit = paste0("Pl:Vx =", nEvents.2, ":", nEvents.1),
                    nPlac       = nEvents.2,
                    nVacc       = nEvents.1,
                    evalTimeCIR = evalTime,
                    varlogFR    = varlogFR,
                    VE_CIR      = VE,
                    VE_CIR_loCI = VE_loCI,
                    VE_CIR_upCI = VE_upCI,
                    stopCrit_CIR = evalStopCrit( cbind(VE_loCI, VE_upCI), 
                                       lowerVE=lowerVE, upperVE=upperVE )
                 )

      ##  <<< This should no longer be needed. Code was added  >>>  ##
      ##  <<< indide of cumIncRatio to handle this case        >>>  ##
      ##  <<<------------------------------------------------- >>>  ##
      ## If there were any test times at which we had no vaccinee infections then
      ## our variance estimates will not be usable - they will (or should) be NA.
      ## Set variances to NA (for cox model) and set 'stopCrit_xxx' appropriately.
      #ZeroVacc <- ( CIRobj$nVacc == 0 )
      #if ( any(ZeroVacc) ) {
      #    CIRobj$stopCrit_CIR[ ZeroVacc ] <- 
      #        ifelse(boundType=="highEff", TRUE, FALSE)
      #}  

      ## subset to retain only the var.s created above
      CIRobj <- CIRobj[c("infectTotal", "infectSplit", "nPlac", "nVacc", 
                       "evalTimeCIR", "VE_CIR", "VE_CIR_loCI", "VE_CIR_upCI",
                       "varlogFR", "stopCrit_CIR") ]

      ## lagged CIR, if requested
      if (laggedMonitoring) {
          cumInc_lagOut <- 
              cumIncRatio( censDatList_lag, times="last", 
                           alphaLevel = alphaLevel,
                           randFraction = randFraction)

          ## subset/reorder columns, then rename (don't need 'futime' here, it'll be
          lagCIRobj <- 
              transform( cumInc_lagOut,
                  infectTotal_lag = nEvents,
                  infectSplit_lag = paste0("Pl:Vx =", nEvents.2, ":", nEvents.1),
                  nPlac_lag       = nEvents.2,
                  nVacc_lag       = nEvents.1,
                  varlogFR_lag    = varlogFR,
                  VE_lagCIR       = VE,
                  VE_lagCIR_loCI  = VE_loCI,
                  VE_lagCIR_upCI  = VE_upCI,
                  stopCrit_lagCIR = evalStopCrit( cbind(VE_loCI, VE_upCI),
                                        lowerVE=lowerVE, upperVE=upperVE )
                  )

          ## if zero infections in either group, then set
          ## 'stopCrit_xxx' appropriately.
          #lagZeroVacc <- ( lagCIRobj$nVacc == 0 )
          #lagZeroPlac <- ( lagCIRobj$nPlac == 0 )

          #if ( any(lagZeroVacc | lagZeroPlac ) ) {
              ## set CI bounds to NA (they could be anything right now)
          #    lagCIRobj$VE_lagCIR_loCI[ lagZeroVacc | lagZeroPlac ] <- NA
          #    lagCIRobj$VE_lagCIR_upCI[ lagZeroVacc | lagZeroPlac ] <- NA

          #    lagCIRobj$stopCrit_lagCIR[ lagZeroVacc ]  <-
          #            ifelse(boundType=="highEff", TRUE, FALSE)

          #    lagCIRobj$stopCrit_lagCIR[ lagZeroPlac ] <-
          #            ifelse(boundType=="nonEff",  TRUE, FALSE)
          #}

          ## subset to retain only the var.s created above
          lagCIRobj <- lagCIRobj[ c("infectTotal_lag", "infectSplit_lag", 
                           "nPlac_lag", "nVacc_lag",
                           "VE_lagCIR", "VE_lagCIR_loCI", "VE_lagCIR_upCI",
                           "varlogFR_lag", "stopCrit_lagCIR") ]
      }
  }


  ## coxph-based hazard ratio, if requested  (same crap, different estimator)
  if (cox) {
      coxHRout <- coxHR( censDatList, alphaLevel = alphaLevel,
                         randFraction = randFraction)

      CoxObj <- transform( coxHRout,
                    infectTotal = nEvents,
                    infectSplit = paste0("Pl:Vx =", nEvents.2, ":", nEvents.1),
                    nPlac       = nEvents.2,
                    nVacc       = nEvents.1,
                    VE_Cox      = VE,
                    VE_Cox_loCI = VE_loCI,
                    VE_Cox_upCI = VE_upCI,
                    stopCrit_Cox = evalStopCrit( cbind(VE_loCI, VE_upCI),
                                       lowerVE=lowerVE, upperVE=upperVE )
                 )

      ##  <<< This should no longer be needed. Code was added  >>>  ##
      ##  <<< indide of cumIncRatio to handle this case        >>>  ##
      ##  <<<------------------------------------------------- >>>  ##
      ## If there were any test times at which we had no vaccinee infections then
      ## our variance estimates will not be usable - they will (or should) be NA.
      ## Set variances to NA (for cox model) and set 'stopCrit_xxx' appropriately.
      #ZeroVacc <- ( CoxObj$nVacc == 0 )
      #if ( any(ZeroVacc) ) {
      #    CoxObj$stopCrit_Cox[ ZeroVacc ] <-
      #            ifelse(boundType=="highEff", TRUE, FALSE)
          ## set CI bounds to NA (they could be anything right now)
      #    CoxObj$VE_Cox_loCI[ ZeroVacc ] <- NA
      #    CoxObj$VE_Cox_upCI[ ZeroVacc ] <- NA
      #}

      ## subset to retain only the var.s created above
      CoxObj <- CoxObj[c("infectTotal", "infectSplit", "nPlac", "nVacc",
                       "VE_Cox", "VE_Cox_loCI", "VE_Cox_upCI", "stopCrit_Cox") ]

      if (laggedMonitoring) {
          coxHR_lagOut <- coxHR( censDatList_lag, alphaLevel = alphaLevel,
                                 randFraction = randFraction)

                lagCoxObj <- transform( coxHR_lagOut,
                    infectTotal_lag = nEvents,
                    infectSplit_lag = paste0("Pl:Vx =", nEvents.2, ":", nEvents.1),
                    nPlac_lag       = nEvents.2,
                    nVacc_lag       = nEvents.1,
                    VE_lagCox       = VE,
                    VE_lagCox_loCI  = VE_loCI,
                    VE_lagCox_upCI  = VE_upCI,
                    stopCrit_lagCox = evalStopCrit( cbind(VE_loCI, VE_upCI),
                                          lowerVE=lowerVE, upperVE=upperVE )
                 )

          ## if zero infections in either group, then set
          ## 'stopCrit_xxx' appropriately.
          #lagZeroVacc <- ( lagCoxObj$nVacc == 0 )
          #lagZeroPlac <- ( lagCoxObj$nPlac == 0 )

          #if ( any(lagZeroVacc | lagZeroPlac ) ) {
          #    ## set CI bounds to NA (they could be anything right now)
          #    lagCoxObj$VE_lagCox_loCI[ lagZeroVacc | lagZeroPlac ] <- NA
          #    lagCoxObj$VE_lagCox_upCI[ lagZeroVacc | lagZeroPlac ] <- NA

          #    lagCoxObj$stopCrit_lagCox[ lagZeroVacc ]  <- 
          #            ifelse(boundType=="highEff", TRUE, FALSE)

          #    lagCoxObj$stopCrit_lagCox[ lagZeroPlac ] <- 
          #            ifelse(boundType=="nonEff",  TRUE, FALSE)
          #}


          ## subset to retain only the var.s created above
          lagCoxObj <- lagCoxObj[,c("infectTotal_lag", "infectSplit_lag", 
                           "nPlac_lag", "nVacc_lag",
                           "VE_lagCox", "VE_lagCox_loCI", "VE_lagCox_upCI",
                           "stopCrit_lagCox") ]
      }
  }

  ## If both cir and cox were specified, then remove columns that are in common.
  ## Do the same for lagged versions, if needed
  if (cir && cox) {
      CoxObj <- CoxObj[ !( names(CoxObj) %in% names(CIRobj) ) ]

      if (laggedMonitoring) {
          lagCoxObj <- lagCoxObj[ (! names(lagCoxObj) %in% names(lagCIRobj) )]
      }
  }


  ##  *** now combine results across estimands  *** 

  ## names of output objects from each estimator
  objNames <- c("CIRobj", "CoxObj", "lagCIRobj", "lagCoxObj")

  ## names of the objects that were computed/created
  objsComputed <- objNames[ indEst ]

  ## It's harder to cbind data.frames that vectors, because they use different
  ## underlying approaches.  The issue is if you're cbind()ing objects that are
  ## conditionally included like this:
  ##    cbind( if (indicator1) a, if (indicator2) b, if (indicator3) c) 
  ##   if a, b, c are vectors it works fine,  if any are data.frames it doesn't
  ##   (unless all indicators are TRUE )
  ##
  ## Anyhoo, that's why I'm using the parse/eval approach to do this
  ## There are other, uglier ways too.  I've not thought of anything cleaner yet 

  cbind_call <- paste0("cbind(", paste(objsComputed, collapse=","), ")")
  summObj <- eval( parse( text = cbind_call ) )

  ## add 'test' and 'testTimes' onto the data.frame
  summObj <- cbind( test = 1:nrow(summObj), testTimes = testTimes, summObj)


  ##  ***** determine if stopping criteria were met at any timepoint  *****
  ##  (overall stopping requires that all individual criteria be met)
  ## identify columns containing stopping criteria evaluation
  w.cols <- which( substr(names(summObj), 1, 8 ) == "stopCrit" ) 
  stopCriteria <- apply( summObj[, w.cols, drop=FALSE], 1, all )


  ## if the stopping criteria was satisfied, subset the output to include
  ## only tests through the *first* time the criteria were met.
  if ( any(stopCriteria) ) {
    boundWasHit <- TRUE
    boundHit <- boundLabel

    stopIndx <- which( stopCriteria )[1]
    summObj  <- summObj[1:stopIndx, ]
    stopTime <- testTimes[ stopIndx ]
    stopInfectCnt <- summObj$infectTotal[ stopIndx ] 
  } else {
    boundWasHit <- FALSE
    boundHit <- NA
    stopTime <- NA
    stopInfectCnt <- NA
  }

  return(
      list( boundWasHit   = boundWasHit,
            boundHit      = boundHit,
            boundType     = boundType, 
            stopTime      = stopTime, 
            stopInfectCnt = stopInfectCnt,
            summObj       = summObj ) 
        )
}
