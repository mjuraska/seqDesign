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
                    ## aveVE should not be specified if 'veByPeriod' is specified
                    aveVE=NULL,
                    ## option 'manual' added to use when using 'veByPeriod', though in the
                    ## code I'll probably just ignorethis variable when veByPeriod is set
                    VEmodel=c("half", "constant", "manual"),
                    vePeriods,
                    ## 'veByPeriod' is used to specify VEs for each period in 'vePeriods', for
                    ## each treatment group.  It can be specified in two ways, one is by providing
                    ## values for 'fullVE' and then providing relative values for each vePeriod
                    ## e.g. if length(vePeriods)=4 (so 3 periods)
                    ##   veByPeriod=list(fullVE=c(0,.6),
                    ##                   C1 = c(1, 1, 1),
                    ##                   T1 = c(.5, .75, 1) )
                    ## Or by directly specifying the VEs without using fullVE, 
                    ## e.g. 
                    ##   veByPeriod=list(C1 = c(0, 0, 0),
                    ##                   T1 = c(.3, .45, 0.6) )
                    ## The first component (after 'fullVE' - if it exists) must be the 
                    ## control/placebo arm, followed by the treatment arms.  The control arm should
                    ## be named 'C1' and the treatment arms 'T1', 'T2', ...
                    ## If fullVE is specified, is should have the same length as the number of arms
                    veByPeriod = NULL,
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

if ( !is.null(veByPeriod) ) {
  ## make a copy then, do some checks
  veBP <- veByPeriod

  if ( !is.null(aveVE) )
    stop("Arguments 'veByPeriod' and 'aveVE' cannot both be specified\n\n")

  if (!is.null(veByPeriod$fullVE) ) {
    fullVE <- veByPeriod$fullVE

    ## remove this component
    veBP$fullVE <- NULL
  
    if ( length(fullVE) != length(veBP) )
      stop("The length of veByPeriod$fullVE must match the number of other components\n",
           "specified in list veByPeriod\n\n")

    ## multiply fullVE by the values provided
    veBP[] <- lapply(1:length(fullVE), function(i, fVE, vbp) {
                     ## multiply the i-th fullVE value by the relative VEs from i-th arm
                     fVE[i]*vbp[[i]] }, fVE=fullVE, vbp=veBP )

  } else {
    ## 'fullVE' was not specified, derive it as the max of the values from each arms
    fullVE <- sapply(veBP, max) 
  }

}

## if 
if (is.null( veByPeriod) ) {
  VEmodel <- match.arg(VEmodel)
} else {
  if (!missing(VEmodel) && VEmodel != "manual")
    cat("Note: You have specified 'veByPeriod' but not set argument 'VEmodel'\n",
        "to 'manual'. The value of 'VEmodel' you provided will be ignored\n\n")
  VEmodel <- "manual" 
}
  

## total number of trial participants
Nppt = sum(N)

if ( !is.null(veByPeriod) ) {
  if (! "C1" %in% names(veBP) )
    stop("There must be a control arm named 'C1' specified in'veByPeriod'\n\n")
 
  ## number of vaccine arms
  nVaccArms = length(veBP) - 1

  if ( any(names(veBP) != c("C1", paste0("T",1:nVaccArms)) ) ) 
    stop("The components of veByPeriod must be named 'C1', 'T1', 'T2', ...\n\n")

  ## the null VE will be the first VE of "C1" 
  nullVE <- veBP[["C1"]][1]
} else {
  ## verify whether length of 'N' = length of 'aveVE'
  if ( length(aveVE) != length(N) ) 
     stop( "Length of 'aveVE' does not match the number of treatment arms given in 'N'.\n" )

  ## VE for placebo arm
  nullVE = aveVE[1]

  ## VE for vaccine arms
  aveVE = aveVE[-1] 

  ## number of vaccine arms
  nVaccArms = length(aveVE)
}

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

if( VEmodel!="constant" && is.null(veByPeriod)) {
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
  
} else if (VEmodel=="constant") {
  vaccEff <- aveVE
}

if (is.null(veByPeriod) ) {
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
} else {
  VEs <- veBP[[-1]]
}

## 'infecRateTbl' contains information on relative infection rates (hazard 
##  ratios) for each treatment.  Please use "Inf" rather than NA to represent 
##  intervals that continue indefinitely.

if (is.null(veByPeriod) ) {
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
} else {

  infecRateList <- lapply(1:length(veBP), function(i,veBP,veP) {
                            data.frame( trt = names(veBP)[i],
                                        start = veP,
                                        end = c(veP[-1], Inf),
                                        relRate= 1 - veBP[[i]] )
                            }, veBP=veBP, veP=vePeriods )
  infecRateTbl  <- do.call(rbind, infecRateList)

}

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
                  fullVE = ifelse(exists("fullVE"), fullVE, vaccEff),
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
      if ( is.null(saveFile) ) {
        infecRate.fmt <- format(infecRate, digits=3, nsmall=3)
        saveFile <- paste0(
           "simTrial_nPlac=", N[1], "_nVacc=", paste(N[-1], collapse="_"), 
           ifelse( exists("aveVE"), 
             paste0("_aveVE=",  paste( round(aveVE,2), collapse="_")),
             paste0("_fullVE=", paste( round(fullVE,2), collapse="_")) ),
           "_infRate=", format(infecRate, digits=3, nsmall=3), ".RData" )

      }
    #save(trialObj, file=file.path(saveDir, saveFile), compress="xz")
    save(trialObj, file=file.path(saveDir, saveFile))

    if (verbose){ 
        cat("Output saved in:\n", file.path(saveDir, saveFile), "\n\n") 
    }
  } 
  return( invisible( trialObj ) )
}

########################### End of simTrial function ###########################



