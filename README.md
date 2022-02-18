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

Simple examples:
```{r}
library(aLDG)

#===== bivariate dependence measure ======#
x = rnorm(100)
# independence
bidep(x, rnorm(length(x),0,0.1) ,methods=c('Pearson','TauStar','HSIC','dCor','aLDG')) 
# linear
bidep(x, x + rnorm(length(x),0,0.1) ,methods=c('Pearson','TauStar','HSIC','dCor','aLDG')) 
# nonlinear
ansquad = bidep(x, x^2 + rnorm(length(x),0.1) ,methods=c('Pearson','TauStar','HSIC','dCor','aLDG')) 
# monotone
bidep(x, x^3 + rnorm(length(x),0,0.1) ,methods=c('Pearson','TauStar','HSIC','dCor','aLDG')) 
# nonmonotone
bidep(x, sin(x*4) + rnorm(length(x),0,0.1) ,methods=c('Pearson','TauStar','HSIC','dCor','aLDG')) 

#===== pairwise dependence matrix for multivariate data ======#
# multivariate independent normal
dat = matrix(rnorm(1000),5,200)
matdep(dat,methods=c('Pearson','TauStar','HSIC','dCor','aLDG'), ncores=NULL)
# multivariate dependent normal
A = diag(5)
A[1,2]=A[2,1]=0.7
matdep(A%*%dat,methods=c('Pearson','TauStar','HSIC','dCor','aLDG'), ncores=NULL)
```

## Reference:
Check out our paper for aLDG here:

From local to global gene co-expression estimation using single-cell RNA-seq data.(to be submitted)
