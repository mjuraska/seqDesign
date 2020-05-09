applyStopRules <- 
    function(d, infectionTotals=NULL, testTimes=NULL, 
             boundType = c("nonEff","eff","highEff"), boundLabel=boundType,
             lowerVE=NULL,  upperVE=NULL, alphaLevel=0.05,
             estimand = c("cox","cuminc","combined"),
             lagTime=NULL, returnAll = FALSE, randFraction)
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
  ## infectionTotals - a vector specifying the total number of infections at which 
  ##                   analyses should take place.  The totals (or "counts") are 
  ##                   assumed to be on the lagged scale (i.e. only lagged events
  ##                   counted) if lagTime is specified
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
  ## lagTime     - The amount of "lag-time" before lagged-monitoring begins.
  ##               Only required if lagged monitoring is desired
  ##
  ## estimand    - a character string specifying the estimand(s) to be used in
  ##                 monitoring (can be one of "combined", "cox", and "cuminc")
  ##
  ## returnAll   - Should the results from all test times be returned?  If FALSE
  ##               (the default) only results through the first test at which the
  ##               stopCriteria are met are returned. If TRUE, then results from
  ##               all tests are returned.  This allows for combining results
  ##               from different cohorts outside this function.
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
  ## We need to do this in a way that takes into account whether lagged-monitoring
  ## was specified.  If so, the counts are assumed to be on the lagged scale.
  if ( is.null(testTimes) ) {
    if ( !is.null(infectionTotals) ) {
        testTimes <- getInfectionTimes(d, cnts=infectionTotals, lagTime=lagTime)  
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
  indEst <- c(cir = cir,  cox = cox )

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

  ## Censor data to apply lagTime, if non-null and and non-zero
  if ( !is.null(lagTime) && length(lagTime)==1 && lagTime>0.001) {
      censDatList <- censorTrial(censDatList, times=lagTime,
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

      ## subset to retain only the var.s created above
      CIRobj <- CIRobj[c("infectTotal", "infectSplit", "nPlac", "nVacc", 
                       "evalTimeCIR", "VE_CIR", "VE_CIR_loCI", "VE_CIR_upCI",
                       "varlogFR", "stopCrit_CIR") ]

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
  }

  ## If both cir and cox were specified, then remove columns that are in common.
  ## Do the same for lagged versions, if needed
  if (cir && cox) {
      CoxObj <- CoxObj[ !( names(CoxObj) %in% names(CIRobj) ) ]
  }


  ##  *** now combine results across estimands  *** 

  ## names of output objects from each estimator
  objNames <- c("CIRobj", "CoxObj")

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
    if ( !returnAll ) {
      stopIndx <- which( stopCriteria )[1]
      summObj  <- summObj[1:stopIndx, ]
      stopTime <- testTimes[ stopIndx ]
      stopInfectCnt <- summObj$infectTotal[ stopIndx ] 
    } else {
      ## when returnAll is TRUE
      stopIndx <- which( stopCriteria )
      stopTime <- testTimes[ stopIndx ]
      stopInfectCnt <- summObj$infectTotal[ stopIndx ] 
    }
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
