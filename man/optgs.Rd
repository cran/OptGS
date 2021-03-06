\name{optgs}
\alias{optgs}
\title{
Finding optimal and balanced group-sequential designs
}
\description{
optgs is used to find a one-sided multi-stage design that balances four optimality criteria for a RCT with normally distributed outcomes
}
\usage{
optgs(delta0 = 0, delta1 = 1/3, J = 2, sigma = 1, sd.known = TRUE,
    alpha = 0.05, power = 0.9, weights = c(0.95, 0, 0, 0.05),
    initial = NULL)
} 
\arguments{
\item{delta0}{mean difference in treatment effect under the null hypothesis (default: 0)}
\item{delta1}{clinically relevant difference used to power the trial (default: 1/3)}
\item{J}{number of stages in the trial (default: 2)}
\item{sigma}{assumed standard deviation of treatment responses (default: 1)}
\item{sd.known}{logical value indicating if sigma will be treated as known; if FALSE, a quantile substitution method will be used to modify the stopping boundaries (default TRUE)}
\item{alpha}{one-sided type-I error rate required (default: 0.05)}
\item{power}{power required (default: 0.9)}
\item{weights}{vector of length 4 giving the weights put on the four optimality criteria (default: c(0.95,0,0,0.05)). See details for more information}
\item{initial}{starting values for the Nelder-Mead algorithm if the user wishes to override the default (default: NULL). Initial values must be specified as a two-dimensional vector where both entries are between -0.5 and 0.5.}
}
\details{optgs uses the extended power-family of group-sequential tests, and searches for the values of the futility and efficacy shape parameters that optimise the specified weighting. A description of the extended power-family and optgs is provided in Wason (2012). The `weights' argument corresponds to the weight put on: 1) the expected sample size at delta=delta0; 2) the expected sample size at delta=delta1; 3) the maximum expected sample size; 4) the maximum sample size (i.e. J*groupsize).}
\value{
\item{groupsize}{the number of patients required per arm, per stage}
\item{futility}{the futility boundaries for the design}
\item{efficacy}{the efficacy boundaries for the design}
\item{ess}{the expected sample size at the delta0; the expected sample size at the delta1; and the maximum expected sample size}
\item{typeIerror}{the actual type-I error rate of the design}
\item{power}{the actual power of the design}
}
\references{Wason, J.M.S. OptGS: an R package for finding near-optimal
  group-sequential designs. Journal of Statistical Software, 66(2),
  1-13. http://www.jstatsoft.org/v66/i02/}

\examples{##Find a three-stage design that minimises the maximum expected sample size.
threestagedeltaminimax=optgs(J=3,weights=c(0,0,1,0)) 
plot(threestagedeltaminimax)
}

