#' Tangent-likelihood function.
#'
#' The function calculates tangent-likelihood given the value of density.
#'
#' @param f values of density function.
#' @param t tangent point. It is positive and close to 0. If 0, the tangent-likelihood is equivalent to log-likelihood.
#' @param p Taylor expansion order, the value can be set up to 3.
#' @return tangent likelihood
#' @noRd
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
#' @noRd

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
#' @return Returns estimates from MTE method.
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
#' @param beta.ini initial estimates of beta. If not specified, LADLasso estimates from \code{rq.lasso.fit()} in \code{rqPen}
#' is used. Otherwise, robust estimators are strongly recommended.
#' @param p Taylor expansion order.
#' @param lambda regularization parameter for LASSO, but not necessary if "adaptive=TRUE".
#' @param adaptive logic argument to indicate if Adaptive-Lasso is used. Default is TRUE.
#' @param t the tangent point. You may specify a sequence of values, so that the function automatically select the optimal one.
#' @param intercept logical input that indicates if intercept needs to be estimated. Default is FALSE.
#' @param penalty.factor can be used to force nonzero coefficients. Default is rep(1, ncol(X)) as in glmnet.
#' @param ... other arguments that are used in \code{glmnet}.
#'
#' @return It returns a sparse vector of estimates of linear regression. It has two types of penalty, LASSO and AdaLasso.
#' Coordinate descent algorithm is used for iteratively updating coefficients.
#' \item{beta}{sparse regression coefficient}
#' \item{fitted}{predicted response}
#' \item{t}{optimal tangent point}
#'
#' @export
#'
#' @examples
#' set.seed(2017)
#' n=200; d=500
#' X=matrix(rnorm(n*d), nrow=n, ncol=d)
#' beta=c(rep(2,6), rep(0, d-6))
#' y=X%*%beta+c(rnorm(150), rnorm(30,10,10), rnorm(20,0,100))
#' output.MTELasso=MTElasso(X, y, p=2, t=0.01)
#' beta.est=output.MTELasso$beta
#'
#' @import glmnet
#' @import rqPen
MTElasso<- function(X, y, beta.ini, p=2, lambda=NULL, adaptive=TRUE, t=0.01, intercept=TRUE, penalty.factor=rep(1,ncol(X)),...){


  if(missing(beta.ini)){
    ### did not pass check due to :::
    # if(intercept==TRUE){
    #   rqfit <- rqPen:::rq.lasso.fit(x=X[,-1], y=as.vector(y), lambda = 0.1, intercept = intercept)
    # }else{
    #   rqfit <- rqPen:::rq.lasso.fit(x=X, y=as.vector(y), lambda = 0.1, intercept = intercept)
    # }
    # beta.ini <- rqfit$coefficients

    rqfit <- rq.pen(x=X, y=as.vector(y), lambda = c(0.2, 0.1))
    if(intercept==TRUE){
      beta.ini <- rqfit$models[[1]]$coefficients[,2]
    }else{
      beta.ini <- rqfit$models[[1]]$coefficients[-1,2]
    }

  }


  if(intercept==TRUE){
    penalty.factor <- c(0, penalty.factor)
    X <- cbind(1, X)
  }else{
    penalty.factor <- penalty.factor
  }

  n<- dim(X)[1]
  d<- dim(X)[2]

  if(missing(p)) p=1

  beta1<- beta.ini
  X1<- X

  delta<- 2
  step<- 0

  id1<- seq(1,ncol(X))


  while((delta>=1e-4) && (step<200) && (length(beta1)>1)){

    step<- step + 1
    #cat("iter",step,"\n")

    beta0<- beta1
    X0<- X1

    d<- length(beta0)

    ## pdf of error
    e<- y-X0%*%beta0
    u<- dnorm(e, mean=0, sd=1)  # assume normal errors

    ## penalty weight - lambda
    if(adaptive==FALSE & is.null(lambda)){
      stop("Error: If adaptive=FALSE, lambda must be provided!")
    }

    if(adaptive==FALSE){
      bic.lambda=lambda*penalty.factor
    }else{
      if(is.null(lambda)){
        bic.lambda<- log(n)/(n*(abs(beta0))+1e-7)*penalty.factor
      }else{
        #bic.lambda<- lambda*log(n)/(n*(abs(beta0))+1e-7)
        bic.lambda<- lambda/(abs(beta0)+1e-7)*penalty.factor
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
        ind1<- which(u0>=t)
        ind2<- which(u0<t)
        L[ind1]<- log(u0[ind1])
        L[ind2]<- log(t) + 1/t*(u0[ind2]-t)*(p==1 | p==2 | p==3) -1/(2*t^2)*(u0[ind2]-t)^2*(p==2 | p==3) + 1/(3*t^3)*(u0[ind2]-t)^3*(p==3)

        #penalty <- t(as.matrix(bic.lambda))%*%abs(theta00)
        penalty <- sum(bic.lambda*abs(theta00))
        return(-sum(L)+n*penalty)
      }

      # update the estimates
      optim.sol<- optimize(PLt, lower=beta0[j]-2, upper=beta0[j]+2)
      theta00[j]<- optim.sol$min
      optim.PLt<- optim.sol$obj

      #cat("beta",j,"\n")
    }## end of coordinate descent

    beta11<- theta00

    id<- which(abs(beta11)>1e-5)
    id1<- id1[id]
    beta1<- beta11[id]
    X1<- X0[,id]

    penalty.factor <- penalty.factor[id]

    delta<- max(abs(beta11-beta0))
    #delta<- sum((beta11-beta0)^2)

  }# end of while loop
  beta<- rep(0,ncol(X))
  beta[id1]<- beta1

  if(intercept==TRUE){
    names(beta) <- c("Intercept", paste("X", 1:(length(beta)-1), sep = ""))
  }else{
    names(beta) <- paste("X", 1:length(beta), sep = "")
  }
  result<- list(beta, X%*%beta, t)
  names(result)<- c("beta", "fitted", "t")
  return(result)

}# End of function
###########################################################









