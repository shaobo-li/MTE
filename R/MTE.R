#' Tangent-likelihood function.
#'
#' The function calculates tangent-likelihood given the value of density.
#'
#' @param f values of density function.
#' @param t tangent point. It is positive and close to 0. If 0, the tangent-likelihood is equivalent to log-likelihood.
#' @param p Taylor expansion order, the value can be set up to 3.
#'
#' @export
#'
#' @examples
#' set.seed(2017)
#' x=c(rnorm(80), rnorm(20, 10, 10))
#' f=dnorm(x)
#' y=lnt(f, 0.1, 2)
lnt<- function(f, t, p){
  if(p>3) stop("Error: p can only be integer smaller than 4")
  t.vec<- rep(t, length(f))
  f1<- ifelse(f>=t,f,t)
  lntlik<- log(f1)*(f>=t) + (log(t) + (f-t)/t*(p==1 | p==2 | p==3) - (f-t)^2/(2*t^2)*(p==2 | p==3) + (f-t)^3/(3*t^3)*(p==3))*(f<t)
  return(lntlik)
}
###########################################################

#' Derivatives of Tangent Likelihood w.r.t. beta
#'
#' @param X design matrix
#' @param u density values
#' @param e eesidual
#' @param p Taylor expansion order
#' @param t0 tangent point
#'
#' @return
#' It returns to first and second derivatives of tangent-likelihood w.r.t. beta.
#' This function is called in MTE main function.
#'
#' @keywords internal

der.LtLik<- function(X, u, e, p, t0){
  n<- dim(X)[1]
  d<- dim(X)[2]

  ## pdf of error
  # e<- y-X%*%beta0
  # u<- dnorm(e, mean=0, sd=s0)

  ## 1st derivative of Lt
  dl<- matrix(0,n,d)
  ind1<- which(u>=t0)
  ind2<- which(u<t0)
  dl[ind1,]<- e[ind1]*X[ind1,]
  dl[ind2,]<- (1-(1-u[ind2]/t0)^p)*e[ind2]*X[ind2,]

  ## 2nd derivative of Lt -- Analytical form
  ddl<- matrix(0,d,d)
  for(i in 1:n){
    ddlt<- matrix(0,d,d)
    if(u[i]>=t0){
      ddlt<- -X[i,]%*%t(X[i,])
    }else{
      ddlt<- (1-(1-u[i]/t0)^p-p*(1-u[i]/t0)^(p-1)*u[i]/t0*e[i]^2)* (-X[i,]%*%t(X[i,]))
    }
    ddl<- ddl+ddlt
  }

  return(list(der1=dl, der2=ddl))
}
#################################################

#' Maximum Tangent-likelihood Estimation
#'
#' It estimates linear regression coefficient using MTE.
#' The function produces robust estimates of linear regression. Outliers and contamination would be downweighted.
#' It is robust to Gaussian assumption of the error term. Initial estimates need to be provided.
#'
#' @param y the response vector
#' @param X design matrix
#' @param beta.ini initial value of estimates, could be from OLS.
#' @param t the tangent point. You may specify a sequence of values, so that the function automatically select the optimal one.
#' @param p Taylor expansion order, up to 3.
#' @param intercept logical input that indicates if intercept needs to be estimated. Default is FALSE.
#'
#'
#' @return
#' \item{beta}{the regression coefficient estimates}
#' \item{fitted.value}{predicted response}
#' \item{t}{the optimal tangent point through data-driven method}
#'
#' @export
#'
#'
#' @examples
#' set.seed(2017)
#' n=200; d=4
#' X=matrix(rnorm(n*d), nrow=n, ncol=d)
#' beta=c(1, -1, 2, -2)
#' y=-2+X%*%beta+c(rnorm(150), rnorm(30,10,10), rnorm(20,0,100))
#' beta0=beta.ls=lm(y~X)$coeff
#' beta.MTE=MTE(y,X,beta0,0.1,2, intercept=TRUE)$beta
#' cbind(c(-2,beta), beta.ls, beta.MTE)
#'
#' @import stats
#'
MTE<- function(y, X, beta.ini, t, p, intercept=FALSE){

  if(intercept==TRUE) X=cbind(1, X)
  n<- dim(X)[1]
  d<- dim(X)[2]

  if(missing(beta.ini)){
    stop("Error: must provide initial estimates")
  }

  if(missing(t)) t=seq(0,0.1,0.005) else t=t

  if(missing(p)) p=1 else p=p

  beta1<- beta.ini
  X0<- X

  delta<- 2
  step<- 0

  while((delta>=1e-4) && (step<20) && (length(beta1)>1)){

    step<- step + 1
    #cat("iter",step,"\n")

    beta0<- beta1

    ## pdf of error
    e<- y-X0%*%beta0
    u<- dnorm(e, mean=0, sd=1)  # assume normal errors

    ## choose tuning parameter t
    if(min(u)>1e-4){
      t0=0
    }else{
      if(d<n/5 & length(t)>1){
        det.V<- rep(0,length(t))
        for(k in 1:length(t)){

          der<- der.LtLik(X=X0, u=u, e=e, p=p, t0=t[k])
          dl<- der$der1   # 1st derivative of Lt
          ddl<- der$der2  # 2nd derivative of Lt

          I<- ddl/n

          # Asymptotic variance Sigma2
          Sigma2<- cov(dl)

          # Calculate det(V)
          V<- solve(I+diag(0.0001,d))%*%Sigma2%*%solve(I+diag(0.0001,d))
          det.V[k]<- det(V)

        }
        #plot(t,det.V,type='l')
        t0<- mean(t[which(det.V==min(det.V))])
      }else{
        t0<- mean(t)
      }
    }

    ## coordinate descent
    theta00<- beta0
    # loops for coordinate descent

    for(j in 1:d){

      # Define univariate function: Lt
      LtFtn<- function(thetaj){
        theta00[j]<- thetaj
        r<- y-X0%*%theta00
        u0<- dnorm(r, mean=0, sd=1)  # assume normal errors
        L<- rep(0,length(u0))
        ind1<- which(u0>=t0)
        ind2<- which(u0<t0)
        L[ind1]<- log(u0[ind1])
        L[ind2]<- log(t0) + 1/t0*(u0[ind2]-t0)*(p==1 | p==2 | p==3) -1/(2*t0^2)*(u0[ind2]-t0)^2*(p==2 | p==3) + 1/(3*t0^3)*(u0[ind2]-t0)^3*(p==3)

        return(sum(L))
      }

      # Define univariate function: PLt
      PLt<- function(thetaj){
        theta00[j]<- thetaj
        PLtf<- LtFtn(thetaj)
        return(-PLtf)
      }

      # update the estimates
      optim.sol<- optimize(PLt, lower=beta0[j]-2, upper=beta0[j]+2)
      theta00[j]<- optim.sol$min
      optim.PLt<- optim.sol$obj


    }## end of coordinate descent

    beta1<- theta00

    delta<- max(abs(beta1-beta0))

  }# end of while loop

  result<- list(X%*%beta1, beta1, t0)
  names(result)<- c("fitted.value", "beta", "t")
  return(result)
}
#########################################################

#' MTE-Lasso estimator
#'
#' MTELasso is the penalized MTE for robust estimation and variable selection for linear regression.
#' It can deal with both fixed and high-dimensional settings.
#'
#' @param y response vector.
#' @param X design matrix, standardization is recommended.
#' @param beta.ini initial estimates of beta. Using unpenalized MTE or LAD is recommended under high-dimensional setting.
#' @param p Taylor expansion order.
#' @param lambda regularization parameter for LASSO, but not necessary if "adaptive=TRUE".
#' @param adaptive logic argument to indicate if Adaptive-Lasso is used. Default is TRUE.
#' @param method it can be ("MTE", "MLE"). The default is MTE. If MLE, classical LASSO is used.
#' @param t the tangent point. You may specify a sequence of values, so that the function automatically select the optimal one.
#' @param intercept logical input that indicates if intercept needs to be estimated. Default is FALSE.
#' @param ... other arguments that are used in function "adalasso()" that is called form parcor package.
#'
#' @return It returns a sparse vector of estimates of linear regression. It has two types of penalty, LASSO and AdaLasso.
#' Coordinate descent algorithm is used for interatively updating coefficients.
#' \item{beta}{sparse regression coefficient}
#' \item{fitted}{predicted response}
#' \item{t}{optimal tangent point}
#'
#' @export
#'
#' @examples
#' set.seed(2017)
#' n=200; d=50
#' X=matrix(rnorm(n*d), nrow=n, ncol=d)
#' beta=c(rep(2,6), rep(0, 44))
#' y=X%*%beta+c(rnorm(150), rnorm(30,10,10), rnorm(20,0,100))
#' beta0=MTE(y, X, rep(0,50), 0.1, 2)$beta
#' output.MTELasso=MTElasso(y,X, p=2, beta.ini=beta0, t=seq(0, 0.1, 0.01), method="MTE")
#' beta.est=output.MTELasso$beta
#'
#' @import glmnet
MTElasso<- function(y, X, beta.ini, p, lambda, adaptive=T, t, method="MTE", intercept=FALSE, ...){

  if(intercept==TRUE) X=cbind(1, X)

  if(method=="MTE"){
    if(missing(t)){
      t=exp(seq(-7, -2, length.out = 20))
    }else{
      t=t
    }
  }
  if(method=="L2")  t=2
  if(method=="MLE"){

    return(glmnet(x=X, y=y, ...))

  }else{

    n<- dim(X)[1]
    d<- dim(X)[2]

    if(missing(beta.ini)){
      stop("Error: must provide initial estimates")
    }

    if(missing(p)) p=1

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

      ## pdf of error
      e<- y-X0%*%beta0
      u<- dnorm(e, mean=0, sd=1)  # assume normal errors


      ## choose tuning parameter t
      if(min(u)>1e-4){
        t0=0
      }else{
        if(d<n/5 & length(t)>1){
          det.V<- rep(0,length(t))
          for(k in 1:length(t)){

            der<- der.LtLik(X=X0, u=u, e=e, p=p, t0=t[k])
            dl<- der$der1   # 1st derivative of Lt
            ddl<- der$der2  # 2nd derivative of Lt

            I<- ddl/n

            # Asymptotic variance Sigma2
            Sigma2<- cov(dl)

            # Calculate det(V)
            V<- solve(I+diag(0.0001,d))%*%Sigma2%*%solve(I+diag(0.0001,d))
            det.V[k]<- det(V)

          }
          #plot(t,det.V,type='l')
          t0<- mean(t[which(det.V==min(det.V))])
        }else{
          t0<- mean(t)
        }
      }


      ## penalty weight - lambda
      if(adaptive==F & missing(lambda)){
        stop("Error: If adaptive=FALSE, lambda must be provided!")
      }

      if(adaptive==F){
        bic.lambda=lambda
      }else{
        if(missing(lambda)){
          bic.lambda<- log(n)/(n*(abs(beta0))+1e-7)
        }else{
          #bic.lambda<- lambda*log(n)/(n*(abs(beta0))+1e-7)
          bic.lambda<- lambda/(abs(beta0)+1e-7)
        }
      }


      ## coordinate descent
      theta00<- beta0
      # loops for coordinate descent

      for(j in 1:d){

        # Define univariate function: Lt
        PLt<- function(thetaj){
          theta00[j]<- thetaj
          r<- y-X0%*%theta00
          u0<- dnorm(r, mean=0, sd=1)  # assume normal errors
          L<- rep(0,length(u0))
          ind1<- which(u0>=t0)
          ind2<- which(u0<t0)
          L[ind1]<- log(u0[ind1])
          L[ind2]<- log(t0) + 1/t0*(u0[ind2]-t0)*(p==1 | p==2 | p==3) -1/(2*t0^2)*(u0[ind2]-t0)^2*(p==2 | p==3) + 1/(3*t0^3)*(u0[ind2]-t0)^3*(p==3)

          return(-sum(L)+n*sum(bic.lambda*abs(theta00)))
        }

        # update the estimates
        optim.sol<- optimize(PLt, lower=beta0[j]-2, upper=beta0[j]+2)
        theta00[j]<- optim.sol$min
        optim.PLt<- optim.sol$obj

        #cat("beta",j,"\n")
      }## end of coordinate descent

      beta11<- theta00
      id<- which(abs(beta11)>1e-4)
      id1<- id1[id]
      beta1<- beta11[id]
      #print(beta1)
      X1<- X0[,id]

      delta<- max(abs(beta11-beta0))
      #delta<- sum((beta11-beta0)^2)

      #print(c(theta1,delta))

    }# end of while loop
    beta<- rep(0,ncol(X))
    beta[id1]<- beta1

    result<- list(beta, X%*%beta, t0)
    names(result)<- c("beta", "fitted", "t")
    return(result)
  }



}# End of function
###########################################################









