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
n=200; d=50
X=matrix(rnorm(n*d), nrow=n, ncol=d)
beta=c(rep(2,6), rep(0, 44))
y=X%*%beta+c(rnorm(150), rnorm(30,10,10), rnorm(20,0,100))
beta0=MTE(y, X, rep(0,50), 0.1, 2)$beta
output.MTELasso=MTElasso(y,X, p=2, beta.ini=beta0, t=seq(0, 0.1, 0.01), method="MTE")
beta.est=output.MTELasso$beta
```

## References

Qin, Y., Li, S., Li, Y., & Yu, Y. (2017). Penalized maximum tangent likelihood estimation and robust variable selection. arXiv preprint [https://arxiv.org/pdf/1708.05439.pdf](arXiv:1708.05439).
