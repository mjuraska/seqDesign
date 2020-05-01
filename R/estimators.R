## Code for the CIR and cox estimation functions


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


