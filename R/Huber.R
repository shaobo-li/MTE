#' Huber Loss
#'
#' @param r residual, y-Xbeta
#' @param alpha 1/alpha is the huber tuning parameter delta. Larger alpha results in smaller portion of squared loss.
#'
#' @return it returns huber loss that will be called in Huber estimation.
#'
#' @export
#'
#'
#'
#'
huberloss<- function(r, alpha){
  if(missing(alpha) | alpha<=0) stop("Error: alpha must be a positive number")
  index<- (abs(r)<=1/alpha)*1
  sum(r^2*index-(2/alpha*abs(r)-1/alpha^2)*(index-1))
}
########################

#' Huber estimation for linear regression
#'
#' This function produces Huber estimates for linear regression. Initial estimates is required.
#' Currently, the function does not support automatic selection of huber tuning parameter.
#'
#' @param y the response vector
#' @param X design matrix
#' @param beta.ini initial value of estimates, could be from OLS.
#' @param alpha 1/alpha is the huber tuning parameter delta. Larger alpha results in smaller portion of squared loss.
#' @param intercept logical input that indicates if intercept needs to be estimated. Default is FALSE.
#'
#' @return
#' \item{beta}{the regression coefficient estimates}
#' \item{fitted.value}{predicted response}
#' \item{iter.steps}{iteration steps.}
#'
#' @export
#'
#' @examples
#' set.seed(2017)
#' n=200; d=4
#' X=matrix(rnorm(n*d), nrow=n, ncol=d)
#' beta=c(1, -1, 2, -2)
#' y=-2+X%*%beta+c(rnorm(150), rnorm(30,10,10), rnorm(20,0,100))
#' beta0=beta.ls=lm(y~X)$coeff
#' beta.huber=huber.reg(y, X, beta0, 2, intercept=TRUE)$beta
#' cbind(c(-2,beta), beta.ls, beta.huber)
huber.reg<- function(y, X, beta.ini, alpha, intercept=FALSE){

  if(intercept==TRUE) X=cbind(1, X)
  n<- dim(X)[1]
  d<- dim(X)[2]

  if(missing(beta.ini)){
    stop("Error: must provide initial estimates")
  }

  beta1<- beta.ini
  X0<- X

  delta<- 2
  step<- 0

  while((delta>=1e-4) && (step<20)){

    step<- step + 1
    #cat("iter",step,"\n")

    beta0<- beta1

    e<- y-X0%*%beta0
    loss<- huberloss(e, alpha)

    ## coordinate descent
    for(j in 1:d){

      # Define univariate function: Lt
      uni.fun<- function(betaj){
        beta0[j]<- betaj
        r<- y-X0%*%beta0
        return(huberloss(r=r, alpha = alpha))
      }

      # update the estimates
      beta0[j]<- optimize(uni.fun, lower=beta0[j]-1, upper=beta0[j]+1)$min

    }## end of coordinate descent


    delta<- max(abs(beta1-beta0))
    beta1<- beta0

  }# end of while loop

  result<- list(step, beta1, X%*%beta1)
  names(result)<- c("iter.steps", "beta", "fitted")
  return(result)
}
###############################################

#' Huber-Lasso estimator
#'
#' This function is L1 penalized Huber estimator for linear regression under both fixed and high-dimensional settings.
#' Currently, the function does not support automatic selection of huber tuning parameter.
#'
#' @param y response vector.
#' @param X design matrix, standardization is recommended.
#' @param beta.ini initial estimates of beta. If not specified, LADLasso estimates from \code{rq.lasso.fit()} in \code{rqPen}
#' is used. Otherwise, robust estimators are strongly recommended.
#' @param lambda regularization parameter of Lasso or adaptive Lasso (if adaptive=TRUE).
#' @param alpha 1/alpha is the huber tuning parameter. Larger alpha results in smaller portion of squared loss.
#' @param adaptive logical input that indicates if adaptive Lasso is used. Default is TRUE.
#' @param intercept logical input that indicates if intercept needs to be estimated. Default is FALSE.
#' @param penalty.factor can be used to force nonzero coefficients. Default is rep(1, ncol(X)) as in glmnet.
#'
#' @return
#' \item{beta}{the regression coefficient estimates.}
#' \item{fitted}{predicted response.}
#' \item{iter.steps}{iteration steps.}
#'
#' @export
#'
#' @examples
#' set.seed(2017)
#' n=200; d=50
#' X=matrix(rnorm(n*d), nrow=n, ncol=d)
#' beta=c(rep(2,6), rep(0, 44))
#' y=X%*%beta+c(rnorm(150), rnorm(30,10,10), rnorm(20,0,100))
#' output.HuberLasso=huber.lasso(X,y, adaptive=TRUE)
#' beta.est=output.HuberLasso$beta
#'
huber.lasso<- function(X, y, beta.ini, lambda, alpha =2, adaptive=TRUE, intercept=FALSE, penalty.factor=rep(1,ncol(X))){

  if(intercept==TRUE) X=cbind(1, X)
  n<- dim(X)[1]
  d<- dim(X)[2]

  if(missing(beta.ini)){
    rqfit <- rqPen::rq.lasso.fit(x=X, y=as.vector(y), lambda = 0.1, intercept = intercept)
    beta.ini <- rqfit$coefficients
  }

  beta1<- beta.ini
  X1<- X

  delta<- 2
  step<- 0
  id1<- seq(1,ncol(X))

  while((delta>=1e-4) && (step<20) && (length(beta1)>1)){

    step<- step + 1
    #cat("iter",step,"\n")

    beta0<- beta1
    X0<- X1
    d<- length(beta0)
    e<- y-X0%*%beta0

    ## penalty weight - lambda
    if(adaptive==FALSE & missing(lambda)){
      stop("Error: If adaptive=FALSE, lambda must be provided!")
    }

    if(adaptive==FALSE){
      bic.lambda=lambda*penalty.factor
    }else{
      if(missing(lambda)){
        bic.lambda<- log(n)/(n*(abs(beta0))+1e-7)*penalty.factor
      }else{
        bic.lambda<- lambda/((abs(beta0))+1e-7)*penalty.factor
      }
    }

    ## coordinate descent
    for(j in 1:d){

      # Define univariate function: Lt
      uni.fun<- function(betaj){
        beta0[j]<- betaj
        r<- y-X0%*%beta0
        l1huber<- huberloss(r=r, alpha = alpha) + n*sum(bic.lambda*abs(beta0))
        return(l1huber)
      }

      # update the estimates
      beta0[j]<- optimize(uni.fun, lower=beta0[j]-2, upper=beta0[j]+2)$min

    }## end of coordinate descent

    delta<- max(abs(beta1-beta0))

    beta11<- beta0
    id<- which(abs(beta11)>1e-4)
    id1<- id1[id]
    beta1<- beta11[id]
    X1<- X0[,id]
    penalty.factor <- penalty.factor[id]

    #delta<- sum((beta11-beta0)^2)

    #print(c(theta1,delta))

  }# end of while loop
  beta<- rep(0,ncol(X))
  beta[id1]<- beta1

  result<- list(step, beta, X%*%beta)
  names(result)<- c("iter.steps", "beta", "fitted")
  return(result)

}# End of function






