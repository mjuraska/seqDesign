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



