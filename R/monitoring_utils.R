## Functions used by monitorTrial 


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
## ------------------------------------------------------------
## Functions for continuous monitoring harm boundary:

##  Input:
##    Bounds: Vector of stopping bounds which correspond to the infection totals
##           1,2,3,... For infection totals where stopping is not permitted (for
##           small infection totals, e.g. 1-7) the 'Bound' should be set to NA.
##
##    p: The randomization probability for the active treatment, when only
##       considering it and the placebo group (in cases when other groups exist).
##       If there are 3 active groups and 1 placebo group, and the randomization
##       prob.s are: 2/9, 2/9, 2/9, 3/9, (the last  being for the placebo group),
##       then the null.p value for generating bounds for one active vs. placebo
##       would be:   2/9 /(2/9 + 3/9) = 2/5
##
##  Function pNS creates an object containing the values of the function P(n,s)
##  defined by Breslow (1970, JASA) as, for 0 <= s <= n <= N,  the probability
##  of the binomial random walk S_n reaching S_n = s without "absorption" into
##  the rejection region (i.e. without hitting any of the stopping bounds).
##
##  Notation:
##    N   = max number of samples
##    S_n = sum{ x_i, i=1,..,n }, and the {x_i} are iid bernoulli(p)
##
##  To calculate P(n,s), we need to specify the probability 'p' the value 'N',
##  and a vector 'B' of length N that specifies the stopping boundaries
##  associated with the values of 'n' from 1:N.  Hence, B[1] will be the stopping
##  value for n=1, B[2] the stopping value for n=2, ..., and B[N] the stopping
##  value for n=N.
##
##  The vector of boundaries should contain NAs for the initial components
##  up until you want to allow stopping (e.g. for n=1,2,3,...,8(?)
##  how far you go will depend on the value of 'p' used.
##
##  The object is a list with components:
##    Bounds - data.frame with columns 'n' = 1:N and
##              'StoppingBound'= contents of input argument 'Bound'
##    p - value of 'p' passed to function
##    N - value of 'N' passed to function
##    totalStopProb - the cumulative probability of stopping by the N-th infection
##   'Stop' - numeric vector of length N, with Stop[i] giving the prob.
##            of hitting the stopping bound (for first time) at n=i
##
##   If argument 'returnPns' is TRUE, then 'Pns' is included in the output obj.
##       Pns is not normally needed, which is why it's been made optional.
##   Pns =  list of length N, each sublist containing a numeric vector.
##            Pns[[ i ]] is a numeric vector of length i+1 containing values
##            of P(n,s), for n=i and s in 0,..,i+1
##
##          Note that values of s >= Bound[i] are zero as they are part of
##          the rejection region (hence there is zero probability of reaching
##          them without hitting the boundary).
##
pNS <- function(Bounds, p=.5, returnPns=FALSE)
{
  N <- length(Bounds)

  if ( all(is.na(Bounds)) )
    stop("The vector provided for argument 'Bounds' contains only NAs.\n",
         "Exiting...\n\n")

  if ( any(Bounds > (1:N), na.rm=TRUE) )
      stop("The bounds provided do not appear to correspond to the infection",
           " totals (1,2,...,N).\n", "One of more of the bounds are larger",
           "than their corresponding infection total.  Exiting...\n\n")

  ## create 'Pns' and 'Stop'. Note that Stop has initial values of 0.
  Stop <- numeric(N)
  Pns  <- vector("list", length=N)

  ## 'first.bound' is the point at which stopping is first allowed
  ## We checked above to make sure that not all values of Bounds are NAs
  ## so we shouldn't have any problem with this...
  first.bound <- which( !is.na( Bounds ) )[1]

  if (first.bound == 1)
    stop("Cannot stop at the first test. Why would you want to?\n\n")


  ## Initialize the base level (i==1). 'i' indexes the total number of infections
  Pns[[ 1 ]] <- c("0" = (1-p), "1" = p)

  ## loop over remaining infection totals, building up the prob.s as we go
  for (i in 2:N) {

    ## initialize vector to store prob.s of the i+1 possible outcome of the
    ## binomial(i, p) variable (the values are: 0,...,i)
    Pns[[ i ]] <- structure( numeric(i+1), names= as.character(0:i) )

    ## Only need to compute for up to the smaller of i and Bounds[i]-1
    ## And since Bounds[i] has to be <= i, Bounds[i] is always < i, so long
    ## as Bounds[i] exists (i.e. as long as it's not NA)
    max.S <- ifelse( is.na(Bounds[i]), i, Bounds[i]-1)

    Pns[[i]]["0"] <- (1-p)*Pns[[i-1]]["0"]
    for (s in 1:max.S )
    {
      ## 'S' is character version of 's'. Used for indexing.
      S         <- as.character( s )
      S.minus.1 <- as.character( s-1 )

      if (s < i )
        Pns[[i]][ S ] <- p*Pns[[i-1]][S.minus.1] + (1-p)*Pns[[i-1]][S]
      else
        ## done when s == i (only happens if Bounds[i] is NA)
        Pns[[i]][ S ] <- p*Pns[[i-1]][S.minus.1]
    }

    ## 'Stop' is the prob. that we encounter a stopping bound for the first
    ## time at infection count 'i'.  Only computed if stopping is possible,
    ## else remains at initialized value of 0.
    if ( !is.na(Bounds[i]) )

      ## If this is the first infection count at which we allow stopping, then
      ## the stop value must be computed allowing for people using non-standard
      ## bounds (i.e. ones that that *don't* begin with bound[n] = n).
      if (i == first.bound && Bounds[i] < i ) {

        ## stop prob. is same as in the 'normal' case (below) *plus* the prob.
        ## that we were already at/above the number of "success" needed to stop
        ## at the previous infection total (i.e. i-1) (for which we should have
        ## been given a stopping bound, but weren't).
        Stop[i] <- p*Pns[[ i-1 ]][ as.character(Bounds[i]-1) ] +
                     sum( Pns[[ i-1 ]][as.character(Bounds[i]:(i-1))] )
      } else {
        ## if 'i' is not the the first total at which stopping is allowed, or
        ## if it is, but we have a proper bound (i.e. Bounds[i] = i) then we do
        ## we compute as this:
        Stop[ i ] <- p*Pns[[ i-1 ]][ as.character(Bounds[i]-1) ]
      }
  }
  outObj <- list(
              totalStopProb = sum(Stop),
              Stop = Stop,
              Bounds  = data.frame(n=1:N, StoppingBound=Bounds),
              N = N,
              p = p )

  ## If 'returnPns' was set to TRUE then add Pns to the output object
  if ( isTRUE(returnPns) )
      outObj$Pns <- Pns

  return( outObj )
}
####################### end of function 'pNS'#############################


### THIS USES 'nonConstBounds' framework to do constant-bounds.
####################### end of function 'pNS'#############################
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
    return( uniroot(getCumAlpha, interval = c(0.000001, 0.05), tol=1e-7,
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
  out <- pNS(Bound=bounds$vaccInfecBound, p=null.p)
  
  names(out$Bounds)[ names(out$Bounds)=="StoppingBound" ] <- "Nvacc"
  boundOut <- within(out$Bounds, {
                       Nplac <- n-Nvacc
                       RR <- round( Nvacc*(1-null.p)/(Nplac*null.p), digits=2)
                     })
  boundOut <- cbind( boundOut,
                     stopProb=round(out$Stop,4),
                     cumStopProb=round(cumsum(out$Stop),4),
                     alphaVal = bounds$cutoff )

  harmBounds <- boundOut
  #names(harmBounds)[1:3] <- c("N", "V", "P") 
   
  # Rename columns from 'old.names' to 'new.names'
  old.names <- c("n","Nvacc","Nplac")
  new.names <- c("N","V","P")
  names(harmBounds)[ match(old.names, names(harmBounds)) ] <- new.names

  ## reorder columns
  ord.columns <- c("N","V","P","RR","stopProb","cumStopProb","alphaVal")
  harmBounds <- harmBounds[, ord.columns]

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
                       lagTime=NULL,
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
##   lagTime        specifies the time, in weeks (since enrollemnt), that events must 
##                  occur after in order to be included in the analysis.
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

    ## if lagTime has been specified, then left censor the data of follow-up time
    if ( !is.null(lagTime) && lagTime > 0.01 ) {
      dat <- censorTrial(dat, times=lagTime, timeScale="follow-up", type="left")
    }

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

# equation to be solved
# all arguments are scalars to have a single equation
eqEstHR <- function(hr, boundType, nullHR, alpha, nEvents, p){
  confLim <- log(hr) + ifelse(boundType=="eff", 1, -1) * qnorm(1 - alpha / 2) * 
    sqrt((1 / nEvents) * (2 + p * hr / (1 - p) + (1 - p) / (p * hr)))
  return(confLim - log(nullHR))
}

# get the root of 'eqEstHR' wrt its first argument
# all arguments are scalars to solve a single equation
getEstHR <- function(boundType, nullHR, alpha, nEvents, randFrac){
  
  int <- if (boundType=="eff"){ c(0.05, nullHR) } else { c(nullHR, 5) }
  estHR <- uniroot(eqEstHR, interval=int, boundType=boundType, nullHR=nullHR, alpha=alpha, nEvents=nEvents, p=randFrac)
  return(estHR$root)
}

# Computes the HR (group 1 / group 2) estimates at the boundary
#
# boundType - boundary type (either "eff" or "nonEff")
# nullHR    - the null hypothesis value of the HR parameter
# alpha     - a vector of two-sided nominal significance levels
# nEvents   - a vector of the total number of events at each analysis (components match those of 'nEvents')
# randFrac  - randomization fraction for group 1
estHRbound <- function(boundType=c("eff", "nonEff"), nullHR, alpha, nEvents, randFrac){
  if (length(alpha)!=length(nEvents)){
    stop("The arguments 'alpha' and 'nEvents' must be of equal length.")
  }
  
  boundType <- match.arg(boundType)
  
  estHR <- NULL
  for (i in 1:length(alpha)){
    estHR <- c(estHR, getEstHR(boundType, nullHR, alpha[i], nEvents[i], randFrac))
  }
  
  return(estHR)
}
