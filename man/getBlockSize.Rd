\name{getBlockSize}
\alias{getBlockSize}
\title{
Determine block size for use in blocked randomization
}
\description{
\code{getBlockSize} returns the minimum block size (possibly within a specified range) that is compatible with a trial's overall treatment assignment totals.
}
\usage{
getBlockSize(nvec, range=c(0,Inf))
}
\arguments{
  \item{nvec}{vector specifying the number of participants to be assigned to each treatment group.  The vector should have one component per group, so that its length equals number of groups. The sum of \code{nvec} should equal the total enrollment for the trial.
}
  \item{range}{(Optional) vector of length two giving the lower and upper bounds (respectively) on block sizes that the user wishes to consider.  
}
}
\details{
The ordering of the components of \code{nvec} is not important, so using \code{nvec = c(x,y,z)} will produce the same results as using \code{nvec = c(z,x,y)}.

In block randomization one does not necessarily want the smallest block size, which is the reason for the existance of the \code{range} argument.  For example, a trial with a 1:1 randomization allocation between two groups would have a minimum block size of 2, which most people would consider to be too small.  So a typical usage of \code{getBlockSize} would be to use \code{range} to set a minimum acceptable block size, through use of vector of form \code{c(lowerBound, Inf)}.  A large trial should probably have a block size on the order of 10-20 or larger, depending on factors including the total trial size and speed of enrollment, so setting a minimum is a good idea.
}

\value{
An integer or NA.  If the user does not specify \code{range}, then the function will always return an integer, which is the smallest block size compatible with the specified vector of treatment group sizes.  If the user \emph{has} specified the \code{range}, then the function adds the further constraint that the block size must lie in the closed interval given by \code{range} (i.e., the block size must be greater-than-or-equal-to \code{range[1]} and less-than-or-equal-to \code{range[2]}).  If there are no compatible block sizes that lie in the given interval, then an NA is returned.

Note that the value returned is the \strong{minimum} block size that is compatible, not necessarily the only one. Any other compatible block sizes (if any exist) will be integer multiples of the minimum size.  You can check the feasibility of various integer multiples by seeing if they divide evenly into the total trial size (i.e., into the sum of \code{nvec}).
}
\examples{

getBlockSize(nvec = c(375, 375) ) 
## specify a minimum block size of 10 (no maximum)
getBlockSize(nvec = c(375, 375), range = c(10, Inf) ) 

getBlockSize( nvec = c(30, 510, 390) )
## require a minimum block size of 10 and maximum of 30 
## (not possible with this nvec, so function returns NA)
getBlockSize( nvec = c(30, 510, 390), range = c(10, 30) )
}
