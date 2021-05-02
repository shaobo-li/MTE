
#' Least Absolute Deviance Estimator for Linear Regression
#'
#' @param y reponse vector
#' @param X design matrix
#' @param intercept logical input that indicates if intercept needs to be estimated. Default is FALSE.
#'
#' @return coefficient estimates
#' @export
#'
#' @examples
#' set.seed(1989)
#' n=200; d=4
#' X=matrix(rnorm(n*d), nrow=n, ncol=d)
#' beta=c(1, -1, 2, -2)
#' y=-2+X%*%beta+c(rnorm(150), rnorm(30,10,10), rnorm(20,0,100))
#' beta.ls=lm(y~X)$coeff
#' beta.MTE=LAD(y,X,intercept=TRUE)
#' cbind(c(-2,beta), beta.ls, beta.MTE)
#'
LAD<- function(y, X, intercept=FALSE){
  if(intercept==TRUE) X=cbind(1, X)
  beta0=rep(0, ncol(X))
  delta<- 2
  step<- 0
  while((step<=20)&&(delta>1e-3)){

    step<- step+1
    beta1<- beta0

    #for(j in sample(length(beta1))){
    for(j in 1:length(beta1)){

      LAD<- function(betaj){
        beta1[j]<- betaj
        return(sum(abs(y-X%*%beta1)))
      }

      beta.est<- optimize(LAD,lower=beta1[j]-5, upper=beta1[j]+5, tol=1e-5)$minimum
      #beta.est<- nlm(object, 1)$est
      beta1[j]<- beta.est

    }

    delta<- max(abs(beta1-beta0))
    beta0<- beta1

  }

  return(beta0)

}
###############################################################################

#' LAD-Lasso for Linear Regression
#'
#' @param y reponse vector
#' @param X design matrix, standardization is recommended.
#' @param beta.ini initial estimates of beta. Using unpenalized LAD is recommended under high-dimensional setting.
#' @param lambda regularization parameter of Lasso or adaptive Lasso (if adaptive=TRUE).
#' @param adaptive logical input that indicates if adaptive Lasso is used. Default is TRUE.
#' @param intercept logical input that indicates if intercept needs to be estimated. Default is FALSE.
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
#' output.LADLasso=LADlasso(y,X, beta.ini=LAD(y, X), lambda=0.2, adaptive=TRUE)
#' beta.est=output.LADLasso$beta
#'
#' @import quantreg
LADlasso<- function(y, X, beta.ini, lambda, adaptive=TRUE, intercept=FALSE){

  if(intercept==T) X=cbind(1, X)

  n<- nrow(X)
  beta1<- beta.ini
  X1<- X

  delta<- 2
  step<- 0

  id1<- seq(1,ncol(X))

  while((delta>=1e-4) & (step<20) & (length(beta1)>1)){

    step<- step + 1

    beta0<-beta1
    X0<- X1
    d<- length(beta1)

    ## penalty weight - lambda
    if(adaptive==F & missing(lambda)){
      stop("Error: If adaptive=FALSE, lambda must be provided!")
    }

    if(adaptive==F){
      bic.lambda=rep(lambda, d)
    }else{
      if(missing(lambda)){
        bic.lambda<- log(n)/(n*(abs(beta0))+1e-7)
      }else{
        #bic.lambda<- lambda*log(n)/(n*(abs(beta0))+1e-7)
        bic.lambda<- lambda/(abs(beta0)+1e-7)
      }
    }

    ystar<- c(y,rep(0,d))
    Xstar<- rbind(as.matrix(X0),diag(n*bic.lambda))

    beta1<- rq(ystar~ Xstar -1, tau=0.5)$coeff

    beta11<- beta1
    id<- which(abs(beta11)>1e-4)
    id1<- id1[id]
    beta1<- beta11[id]
    X1<- X0[,id]

    delta<- max(abs(beta11-beta0))

  }# end of while loop
  beta<- rep(0,ncol(X))
  beta[id1]<- beta1

  result<- list(beta, X%*%beta, step)
  names(result)<- c("beta", "fitted", "iter.steps")
  return(result)


}# End of function
###############################################################################
