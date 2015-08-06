
#' @title Kriging model with convexity 1D constraint
#' @param design 1-column matrix of the design of experiments
#' @param response a vector containing the output values given by the real function at the design points
#' @param basis.size a value represents the number of the basis functions (descritization of 1D input set)
#' @param covtype an optimal character string specifying the covariance function to be used ("gauss" and "matern3_2" choice)
#' @param coef.cov a value corresponding to the length theta hyper-parameters of covariance function
#' @param coef.var a value specifying the variance parameter
#' @param nugget an optimal value used as nugget effect to solve the numerical inverse matrix problem
#' @import quadrpog
#' kmConvex1D(design=c(0.2, 0.5, 0.9), response=c(3, -5, 8))
kmConvex1D <- function(design, response, 
                       basis.size = dim(design)[1]+2+10, 
                       covtype = "gauss",
                       coef.cov = 0.5*(max(design)-min(design)),
                       coef.var = 10*var(response),
                       nugget = 1e-8) {
  
  if (!is.matrix(design)) design=matrix(design,ncol=1)
  
  n <- nrow(design) # nb de points ? interpoler
  N <- basis.size
  p <- (N + 2) - n
  u <- seq(0, 1, by = 1/N) # vecteur de discrétisation
  delta <- 1/N
  
  sig <- sqrt(coef.var)
  theta <- coef.cov
  
  
  if(covtype=='gauss'){
    # Noyau gaussien du processus Y
    k <- function(x, xp, sig, theta){
      (sig^2)*exp(-(x-xp)^2/(2*theta^2))
    }
    
    # D?riv?e % x
    kp1 <- function(x, xp, sig, theta){
      -(x-xp)/(theta^2)*k(x,xp, sig, theta)
    }
    
    # D?riv?e % xp
    kp2 <- function(x, xp, sig, theta){
      -kp1(x,xp, sig, theta)
    }
    
    
    # Noyau du processus d?riv?e (D?riv?e % x et xp (ou l'inverse))
    
    kpp <- function(x, xp, sig, theta){
      (1/(theta^2))*k(x,xp, sig, theta)*(1-(x-xp)^2/(theta^2))
    }
    
    # D?riv?e 2 fois % x ou 2 fois % xp
    kpp12 <- function(x,xp, sig, theta){
      -kpp(x,xp, sig, theta)
    }
    
    ## Noyau du processus d?riv?e seconde (2 fois % x et 2 fois % xp)
    
    k2pp <- function(x, xp, sig, theta){
      exp(-(x-xp)^2/(2*theta^2))*((sig^2*(x-xp)^4)/theta^8-(6*sig^2*(x-xp)^2/theta^6)+3*(sig^2)/theta^4)
    }
    
    # D?riv?e 3 fois dont 2 fois % x et une fois % xp
    k3pp1 <- function(x, xp, sig, theta){
      exp(-(x-xp)^2/(2*theta^2))*(sig^2*(x-xp))/(theta^4)*(-3+(x-xp)^2/theta^2)
    }
    
    # D?riv?e 3 fois dont 2 fois % xp et une fois % x
    
    k3pp2 <- function(x, xp, sig, theta){
      exp(-(x-xp)^2/(2*theta^2))*(sig^2*(x-xp))/(theta^4)*(+3-(x-xp)^2/theta^2)
    }
    
    
  }
  else if (covtype=='matern3_2'){
    
  }
  
  else stop('covtype', covtype, 'not supported')
  
  
  
  
  
  
  
  
  
  Gamma <- matrix(data = 0, nrow = N+3, ncol = N+3)
  Gamma[1,1] <- k(0,0, sig, theta)
  Gamma[1,2] <- kp2(0,0, sig, theta)
  Gamma[2,1] <- kp1(0,0, sig, theta)
  Gamma[2,2] <- kpp(0,0, sig, theta)
  for(j in 3 : (N+3)){
    Gamma[1,j] <- kpp12(0, u[j-2], sig, theta)
    Gamma[2,j] <- k3pp2(0,u[j-2], sig, theta)
  }
  for(i in 3 : (N+3)){
    Gamma[i,1] <- kpp12(u[i-2], 0, sig, theta)
    Gamma[i,2] <- k3pp1(u[i-2],0, sig, theta)
  }
  for(i in 3 : (N+3)){
    for(j in 3 : (N+3)){
      Gamma[i, j] = k2pp(u[i-2], u[j-2], sig, theta)
    }
  }
  Gamma <- Gamma + nugget * diag(N+3)
  
  
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
  
  
  
  A <- matrix(data = 0, ncol = N+3, nrow = n)
  for(i in 1 : n){ 
    A[i,1] = 1
    A[i,2] = design[i,1]
    A[i,3] = phi0(design[i,1], N)
    A[i,N+3] = phiN(design[i,1], N)
    for(j in 4 : (N+2)){
      A[i,j] = phii(design[i,1], j-3, N)
    }
  }
  
  
  Amat <- diag(ncol(A))
  Amat[1,1] <- 0
  Amat[2,2] <- 0
  Amat <- rbind(A,Amat)
  
  invGamma1 <- chol(Gamma)
  invGamma <- chol2inv(invGamma1)
  
  zetoil <- solve.QP(invGamma,dvec=rep(0,ncol(A)),Amat=t(Amat),bvec=c(response,rep(0,ncol(A))),meq=n)$solution
  
  return(structure(list(zetoil=zetoil,phi0=phi0,phiN=phiN,phii=phii,Amat=Amat,Gamma=Gamma, A=A, D=D,
                        call=list(design=design, response=response,basis.size=basis.size,covtype=covtype, coef.cov=coef.cov,coef.var=coef.var, nugget=nugget))
                   , class = 'kmConvex1D'))
  
}

#' @title Apply basis functions to given data in design space
#' @param model km* model
#' @param newdata data in design space
#' @method Phi1D kmConvex1D
Phi1D.kmConvex1D <- function(model, newdata){
  N <- model$call$basis.size  
  x <- newdata
  
  v <- matrix(0, nrow = length(x), ncol = N+3)
  v[,1] <- 1
  v[,2] <- x
  v[,3] <- apply(t(x), 2, model$phi0, N = N)
  v[,N+3] <- apply(t(x), 2, model$phiN, N = N)
  for(j in 4 : (N+2)){
    v[,j] = model$phii(x, j-3, N)
  }
  return(v)
}