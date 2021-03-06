% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bspGenerics.R
\name{plot.betaStacyProcess}
\alias{plot.betaStacyProcess}
\title{Plots a Beta-Stacy process over its support}
\usage{
\method{plot}{betaStacyProcess}(x, withConfInt = FALSE,
  withPaths = FALSE, conf.level = 0.95)
}
\arguments{
\item{x}{The bsp object to be plotted}

\item{withConfInt}{If true, will calculate and plot approximate 95\% confidence intervals, default FALSE}

\item{withPaths}{If true, will plot 100 sampled CDFs in gray (must also set withConfInt=T), default FALSE}

\item{conf.level}{Confident level for credible interval, default 0.95}
}
\value{
A ggplot object
}
\description{
Plots a Beta-Stacy process over its support
}
\examples{
plot(bsp(c(1,5),c(0.2,0.4),c(2,1)))


}
