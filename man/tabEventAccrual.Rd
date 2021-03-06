% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eventAccrual.R
\name{tabEventAccrual}
\alias{tabEventAccrual}
\title{Tabulate event accrual over time since first enrollment}
\usage{
tabEventAccrual(
  trialData,
  atEvents = NULL,
  atWeeks = NULL,
  prob = 0.5,
  lagTimeMITT = 0,
  lagTimePP = NULL,
  namePP = "pp1",
  na.ub = 0.2
)
}
\arguments{
\item{trialData}{either a list of data frames from \code{\link{simTrial}} (i.e., component \code{trialData} from the output list) or a character string specifying a path to an \code{.RData} file outputted by \code{\link{simTrial}}}

\item{atEvents}{a numeric vector specifying treatment-pooled event counts for which empirical quantiles of time-to-accrual shall be calculated}

\item{atWeeks}{a numeric vector specifying time points (in weeks) since first enrollment for which empirical quantiles of treatment-pooled event counts shall be calculated}

\item{prob}{a numeric value in \eqn{(0, 1)} specifying the probability at which the empirical quantiles across the simulated trials are computed (default is 0.5)}

\item{lagTimeMITT}{a time point (in weeks). Only events with time-to-event \eqn{\ge} \code{lagTimeMITT} are counted in the MITT column (default is 0).}

\item{lagTimePP}{a time point (in weeks). If specified, only PP events with time-to-event \eqn{\ge} \code{lagTimePP} are counted in the PP column.}

\item{namePP}{a character string specifying the name of the column in each data frame in \code{trialData} which indicates membership in the PP cohort (default is "\code{pp1}")}

\item{na.ub}{a numeric value specifying an upper limit on the fraction of simulated trials that do not reach a given event count in \code{atEvents} to still compute the empirical quantile of time-to-accrual. If the fraction of such trials exceeds \code{na.ub}, \code{NA} will be produced.}
}
\value{
A data frame (with at least two columns) of event totals and associated time periods since first enrollment required to accrue the event totals in the MITT cohort. If an PP cohort is specified (via \code{lagTimePP}), a third column is added.
}
\description{
Tabulates side-by-side different event totals and time periods since first enrollment required to accrue the event totals. The user specifies a vector of either event totals or time points since first enrollment, and \code{tabEventAccrual} completes the table.
}
\details{
All time variables use week as the unit of time.

If the user specifies \code{atEvents}, time periods since first enrollment are computed that are needed to observe \code{atEvents} MITT and \code{atEvents} PP events.

If the user specifies \code{atWeeks}, MITT and PP event totals are computed that are observed by \code{atWeeks} weeks since first enrollment.

The function inputs a large number of simulated trial data, and the computed variables (time periods or event totals) are empirical quantiles at probability \code{prob} of the sample distributions (of time periods or event totals). Medians are computed by default.
}
\examples{
simData <- simTrial(N=c(1000, rep(700, 2)), aveVE=seq(0, 0.4, by=0.2), 
                    VEmodel="half", vePeriods=c(1, 27, 79), enrollPeriod=78, 
                    enrollPartial=13, enrollPartialRelRate=0.5, dropoutRate=0.05, 
                    infecRate=0.04, fuTime=156, 
                    visitSchedule=c(0, (13/3)*(1:4), seq(13*6/3, 156, by=13*2/3)),
                    missVaccProb=c(0,0.05,0.1,0.15), VEcutoffWeek=26, nTrials=5, 
                    blockSize=NULL, stage1=78, randomSeed=300)

## user specifies MITT event totals
tabEventAccrual(simData$trialData, atEvents=seq(10, 100, by=10))

## user specifies MITT and PP event totals
tabEventAccrual(simData$trialData, atEvents=seq(10, 100, by=10), lagTimePP=6)

## user specifies time points since first enrollment
tabEventAccrual(simData$trialData, atWeeks=seq(52, 156, by=8), lagTimePP=6)

}
\seealso{
\code{\link{simTrial}}
}
