\name{all.class}
\alias{plot.OptGS}
\alias{print.OptGS}
\title{Different generic functions for class OptGS}
\description{Generic functions for summarising an object of class OptGS}
\usage{
\method{print}{OptGS}(x,... )
\method{plot}{OptGS}(x,ylim=NULL,...)
}
\arguments{
  \item{x}{An output object of class OptGS}
  \item{ylim}{y limits to be passed to plot}
\item{...}{Additional arguments to be passed.}
}
\details{
print.OptGS gives the group-size, stopping boundaries, and operating characteristics of the design

plot.OptGS produces a plot of the expected sample size as the standardised treatment effect differs
}

\value{Screen or graphics output.}