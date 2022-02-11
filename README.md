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
ansindep = bidep(x, rnorm(length(x),0,0.1) ,methods=c('Pearson','TauStar','HSIC','dCor','aLDG')) # independence
ansindep
anslinear = bidep(x, x + rnorm(length(x),0,0.1) ,methods=c('Pearson','TauStar','HSIC','dCor','aLDG')) # linear
anslinear
ansquad = bidep(x, x^2 + rnorm(length(x),0.1) ,methods=c('Pearson','TauStar','HSIC','dCor','aLDG')) # nonlinear
ansquad

#===== pairwise dependence matrix for multivariate data ======#
dat = matrix(rnorm(1000),5,200)
ans = matdep(dat,methods=c('Pearson','TauStar','HSIC','dCor','aLDG'), ncores=NULL) # multivariate normal
ans
```

## Reference:
Check out our paper for aLDG here:
[From local to global gene co-expression estimation using single-cell RNA-seq data.](to be submitted)
