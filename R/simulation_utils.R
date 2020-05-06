## Functions used by simTrial


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
## lagTime: Specifies a lagTime to use for inclusion of infections.  Only infections
##          occurring after the specified lagTime (measured from enrollment) will be
##          counted.
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

getInfectionTimes <- function( d , cnts, arms=NULL, lagTime=NULL, ..., SIMPLIFY=TRUE)
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

  ## if lagTime has been specified, then left censor the data to that follow-up time
  if ( !is.null(lagTime) && lagTime>0.001 ) {
    d <- censorTrial(d, times=lagTime, timeScale="follow-up", type="left")
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



