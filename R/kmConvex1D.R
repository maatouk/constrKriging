
#' @title Kriging model with convexity 1D constraint
#' @param design 1-column matrix of the design of experiments
#' @param response a vector containing the output values given by the real function at the design points
#' @param basis.size a value represents the number of the basis functions (descritization of 1D input set)
#' @param covtype an optimal character string specifying the covariance function to be used ("gauss" and "matern3_2" choice)
#' @param coef.cov a value corresponding to the length theta hyper-parameters of covariance function
#' @param coef.var a value specifying the variance parameter
#' @param nugget an optimal value used as nugget effect to solve the numerical inverse matrix problem
#' @import MASS
#' @import solve.QP
#' @examples
#' kmConvex1D(design=c(0.1, 0.5, 0.9), response=c(10, 5, 9), coef.cov = 0.3)
#' kmConvex1D(design=c(0.1, 0.5,.7, 0.9), response=c(10, 5,7, 9), coef.cov = 0.3)

kmConvex1D <- function(design, response, 
                       basis.size = dim(design)[1]+2+10, 
                       covtype = "gauss",
                       coef.cov = 0.5*(max(design)-min(design)),
                       coef.var = var(response),
                       nugget = 1e-7*sd(response)) {
  
  if (!is.matrix(design)) design=matrix(design, ncol = 1)
  
  #   if (coef.cov=="LOO") {
  #     object=kmConvex1D(design, response, 
  #                       basis.size, 
  #                       covtype ,
  #                       coef.cov = 0.5*(max(design)-min(design)),
  #                       coef.var,
  #                       nugget)
  #     
  #     theta=coef.cov_LOO(object)
  #     
  #     model = NULL
  #     while(is.null(model)) # auto raise nugget if needed
  #       try(model <- kmConvex1D(design, response, 
  #                       basis.size, 
  #                       covtype ,
  #                       coef.cov = theta,
  #                       coef.var,
  #                       nugget=nugget*10))
  #     
  #     return(model)
  #   }
  
  n <- nrow(design)  # number of design points 
  N <- basis.size    # discretization size
  p <- (N + 3) - n   # degree of freedom
  u <- seq(0, 1, by = 1/N)  # discretization vector 
  delta <- 1/N
  
  sig <- sqrt(coef.var)
  theta <- coef.cov
  
  
  if(covtype=='gauss'){
    ## Gaussian covariance kernel 
    k <- function(x, xp, sig, theta){
      (sig^2)*exp(-(x-xp)^2/(2*theta^2))
    }
    ## first partial derivative
    kp1 <- function(x, xp, sig, theta){
      -(x-xp)/(theta^2)*k(x, xp, sig, theta)
    }
    ## second partial derivative
    kp2 <- function(x, xp, sig, theta){
      -kp1(x, xp, sig, theta)
    }
    ## derivation of the Gaussian kernel with respect to the 1er and second variable
    kpp <- function(x, xp, sig, theta){
      (1/(theta^2))*k(x, xp, sig, theta)*(1-(x-xp)^2/(theta^2))
    }
    ## two times derivative with respect to the first or the second variable
    kpp12 <- function(x,xp, sig, theta){
      -kpp(x,xp, sig, theta)
    }
    ## two times derivative with respect to the first and the second variable
    k2pp <- function(x, xp, sig, theta){
      exp(-(x-xp)^2/(2*theta^2))*((sig^2*(x-xp)^4)/theta^8-(6*sig^2*(x-xp)^2/theta^6)+3*(sig^2)/theta^4)
    }
    ## third times derivative (two times with respect to the first variable and one time to the second one)
    k3pp1 <- function(x, xp, sig, theta){
      exp(-(x-xp)^2/(2*theta^2))*(sig^2*(x-xp))/(theta^4)*(-3+(x-xp)^2/theta^2)
    }
    ## third times derivative (one time with respect to the first variable and two times to the second one)
    k3pp2 <- function(x, xp, sig, theta){
      exp(-(x-xp)^2/(2*theta^2))*(sig^2*(x-xp))/(theta^4)*(+3-(x-xp)^2/theta^2)
    }
  }
  
  
  else stop('covtype', covtype, 'not supported')
  
  
  fctGamma=function(.theta){
    Gamma <- matrix(data = NA, nrow = N+3, ncol = N+3)
    Gamma[1, 1] <- k(0, 0, sig, .theta)
    Gamma[1, 2] <- kp2(0, 0, sig, .theta)
    Gamma[2, 1] <- kp1(0, 0, sig, .theta)
    Gamma[2, 2] <- kpp(0, 0, sig, .theta)
    for(j in 3 : (N+3)){
      Gamma[1, j] <- kpp12(0, u[j-2], sig, .theta)
      Gamma[2, j] <- k3pp2(0, u[j-2], sig, .theta)
    }
    for(i in 3 : (N+3)){
      Gamma[i, 1] <- kpp12(u[i-2], 0, sig, .theta)
      Gamma[i, 2] <- k3pp1(u[i-2], 0, sig, .theta)
    }
    for(i in 3 : (N+3)){
      for(j in 3 : (N+3)){
        Gamma[i, j] = k2pp(u[i-2], u[j-2], sig, .theta)
      }
    }
    Gamma <- Gamma + nugget * diag(N+3)
    return(Gamma)
  }
  Gamma=fctGamma(theta)
  
  ## basis functions (double primitive of the hat functions)
  phi0 <- function(x, N){
    ifelse(x <= -delta, -delta*x/2-delta^2/6, ifelse(x >= -delta & x <= 0, 
                                                     (x+delta)^3/(6*delta)-delta*x/2-delta^2/6,
                                                     ifelse(x >= 0 & x <= delta, (delta-x)^3/(6*delta)+delta*x/2-delta^2/6, delta*x/2-1/6*delta^2)))
  }
  phii <- function(x, i, N){
    ifelse(x <= u[i], 0, ifelse(x >= u[i] & x <= u[i+1], (delta+x-u[i+1])^3/(6*delta),
                                ifelse(x >= u[i+1] & x <= u[i+2], delta*(x-u[i+1])+(delta-x+u[i+1])^3/(6*delta),
                                       delta^2+delta*(x-u[i+2]))))
  }
  phiN <- function(x, N){
    phii(x-delta, N-1)
  }
  
  A <- matrix(data = NA, ncol = N+3, nrow = n)
  A[, 1] = 1
  A[, 2] = design[, 1]
  A[, 3] = phi0(design[, 1], N)
  A[, N+3] = phiN(design[, 1], N)
  for(j in 4 : (N+2)){
    A[, j] = phii(design[, 1], j-3, N)
  }
  
  
  Amat <- diag(ncol(A))
  Amat[1, 1] <- 0
  Amat[2, 2] <- 0
  Amat <- rbind(A, Amat)
  
  invGamma <- chol2inv(chol(Gamma))
  
  zetoil <- solve.QP(invGamma, dvec=rep(0, ncol(A)), Amat=t(Amat), bvec=c(response, rep(0, ncol(A))), meq=n)$solution
  
  return(structure(list(zetoil=zetoil, phi0=phi0, phiN=phiN, phii=phii, Amat=Amat, Gamma=Gamma, A=A, D=D, fctGamma=fctGamma, invGamma=invGamma,
                        call=list(design=design, response=response, basis.size=basis.size, covtype=covtype, coef.cov=coef.cov, coef.var=coef.var, nugget=nugget))
                   , class = 'kmConvex1D'))
  
}

#' @title Apply basis functions to given data in design space
#' @param model km* model
#' @param newdata data in design space
#' @method Phi1D kmConvex1D
Phi1D.kmConvex1D <- function(model, newdata){
  N <- model$call$basis.size  
  x <- newdata
  
  v <- matrix(NA, nrow = length(x), ncol = N+3)
  v[, 1] <- 1
  v[, 2] <- x
  v[, 3] <- model$phi0(x, N = N)
  v[, N+3] <- model$phiN(x, N = N)
  for(j in 4 : (N+2)){
    v[, j] = model$phii(x, j-3, N)
  }
  return(v)
}