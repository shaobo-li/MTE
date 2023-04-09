# MTE: Maximum Tangent Likelihood Estimation

## Overview

The package provides several robust estimation methods for linear regression under both fixed and high dimesional settings. The methods include Maximum Tangent Likelihood Estimator (`MTE` and `MTElasso`) (Qin et al., 2017+), Least Absolute Deviance Estimator (`LAD` and `LADlasso`) and Huber estimator (`huber.reg` and `huber.lasso`).

## Installation

``` r 
devtools::install_github("shaobo-li/MTE")
```

## Example

``` r
library(MTE)
set.seed(2017)
n=200; d=500
X=matrix(rnorm(n*d), nrow=n, ncol=d)
beta=c(rep(2,6), rep(0, d-6))
y=X%*%beta+c(rnorm(150), rnorm(30,10,10), rnorm(20,0,100))
output.MTELasso=MTElasso(X, y, p=2, t=0.01)
beta.est=output.MTELasso$beta
```

## References

Qin, Y., Li, S., Li, Y., & Yu, Y. (2017). Penalized maximum tangent likelihood estimation and robust variable selection. [arXiv:1708.05439](https://arxiv.org/pdf/1708.05439.pdf).
