% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corr.R
\name{matdep}
\alias{matdep}
\title{Function to compute the multivariate dependency matrix}
\usage{
matdep(
  data,
  methods = NULL,
  all = FALSE,
  thred = -Inf,
  wd = 1,
  trial = 0,
  ncores = NULL,
  norm = FALSE,
  abs = FALSE,
  info = "",
  dir = "./dat/"
)
}
\arguments{
\item{data:}{feature by sample matrix}

\item{methods:}{a vector of strings, could be subset of c('pearson','kendall','taustar','dcor','hsic','hoeffd','ssd')}

\item{all:}{logical value, if true then compute for all methods in c('pearson','kendall','taustar','dcor','hsic','hoeffd','ssd')}

\item{thred:}{only the data > thred will be used}

\item{wd:}{the coefficient before the bandwidth h = wd*h_n}

\item{qd:}{the quantile window to compute the adaptive bandwidth h_x = sd(x(qd))}

\item{info:}{integer denote the trial number}

\item{ncores:}{number of cores to use}

\item{norm:}{whether scale the value to [0,1]}

\item{abs:}{whether take absolute value of the value}

\item{extrainfo:}{some setting related information}
}
\value{
a list of matrix
}
\description{
Function to compute the multivariate dependency matrix
}
\examples{
# Load example data
ans = matdep(runif(100,5,20), methods='ssd')
}
