# aLDG

Welcome to ``aLDG``!  ``aLDG`` is an R package for the computing bivariate dependence measure, 
including dozens of canonical methods and state-of art methods.  This packaged was built along with 
a newly proposed bivariate dependence measure called averaged Local Density Gap (aLDG).

## Install from Github
This package can be installed with R package devtools:
```{r}
library("devtools")
devtools::install_github("JINJINT/aLDG")
```

## Quick start:

For a simple example to simulate 100 genes and 50 cells of one cell group with gene co-expression:
```{r}
library(aLDG)

#===== bivariate dependence measure ======#
x = rnorm(100)
anslinear = bidep(x, x + rnorm(100,0,0.1) ,methods=c('pearson','taustar','hsic','dcor','aldg')) # compute for selected methods
ansquad = bidep(x, x^2 + rnorm(100,0,0.1) ,methods=c('pearson','taustar','hsic','dcor','aldg')) # compute for selected methods

#===== pairwise dependence matrix for multivariate data ======#
dat = matrix(rnorm(100),5,20)
ans = matdep(dat,methods=c('pearson','taustar','hsic','dcor','aldg')), ncores=NULL) # compute for selected methods
```

## Reference:
Check out our paper for aLDG here:
[From local to global gene co-expression estimation using single-cell RNA-seq data.](to be submitted)
