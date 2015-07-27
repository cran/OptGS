\name{powerfamily}
\alias{powerfamily}
\title{
Finding extended power-family group-sequential designs
}
\description{
powerfamily is used to find a one-sided extended power-family group-sequential design 
}
\usage{
powerfamily(futility = 0, efficacy = 0, delta0 = 0, delta1 = 1/3,
    J = 2, sigma = 1, sd.known = TRUE, alpha = 0.05, power = 0.9)
}
\arguments{
\item{futility}{shape parameter for futility boundaries (default: 0)}
\item{efficacy}{shape parameter for efficacy boundaries (default: 0)}
\item{delta0}{mean difference in treatment effect under the null hypothesis (default: 0)}
\item{delta1}{clinically relevant difference used to power the trial (default: 1/3)}
\item{J}{number of stages in the trial (default: 2)}
\item{sigma}{assumed standard deviation of treatment responses (default: 1)}
\item{sd.known}{logical value indicating if sigma will be treated as known; if FALSE, a quantile substitution method will be used to modify the stopping boundaries (default TRUE)}
\item{alpha}{one-sided type-I error rate required (default: 0.05)}
\item{power}{power required (default: 0.9)}
}
\details{powerfamily uses the extended power-family of group-sequential tests. A description of the extended power-family is provided in Wason (2012). }
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
\examples{
##Find a three-stage design that has shape parameters -0.5 and 0.5.
threestagedesign=powerfamily(J=3,futility=-0.5,efficacy=0.5) 
plot(threestagedesign)
}
