
#' @title Kriging model with bounds constraints
#' @param design 1-column matrix of the design of experiments
#' @param response a vector containing the output values given by the real function at the design points
#' @param basis.size a value represents the number of the basis functions (descritization of 1D input set)
#' @param covtype an optimal character string specifying the covariance function to be used ("gauss" and "matern3_2" choice)
#' @param coef.cov a value corresponding to the length theta hyper-parameters of covariance function
#' @param coef.var a value specifying the variance parameter
#' @param nugget an optimal value used as nugget effect to solve the numerical inverse matrix problem
#' @param lower lower bound constraint
#' @param upper upper bound constraint

#' @examples 
#' model = kmBounded(design=c(0.1, 0.3, 0.5, 0.9), response=c(7, -8, 9, 15), lower=-10, upper = 18)
kmBounded <- function(design, response, 
                      basis.size = dim(design)[1]+2+10, 
                      covtype = "gauss",
                      coef.cov = 0.5*(max(design)-min(design)),
                      coef.var = var(response),
                      lower = min(response)-(max(response)-min(response))*.1,
                      upper = max(response)+(max(response)-min(response))*.1,
                      nugget = 1e-8) {
  
  if (!is.matrix(design)) design=matrix(design,ncol=1)
  
  n <- nrow(design) # nb de points ? interpoler
  N <- basis.size
  p <- (N + 2) - n
  u <- seq(0, 1, by = 1/N) # vecteur de discrétisation
  
  sig <- sqrt(coef.var)
  theta <- coef.cov
  
  
  if(covtype=='gauss'){
    # Noyau gaussien du processus Y
    k <- function(x, xp, sig, theta){
      (sig^2)*exp(-(x-xp)^2/(2*theta^2))
    }
    
  }
  else if (covtype=='matern3_2'){
    
  }
  
  else stop('covtype', covtype, 'not supported')
  
  
  
  
  
  
  
  phi <- function(x){
    ifelse(x >= -1 & x <= 1, 1-abs(x), 0)
  }
  
  phi0 <- function(x, N){
    phi(x*N)
  }
  
  phiN <- function(x, N){
    phi((x-u[N+1])*N)
  }
  
  phii <- function(x, i, N){
    phi((x - u[i+1])*N)
  }
  
  A <- matrix(data = 0, ncol = N+1, nrow = n)
  for(i in 1 : n){ 
    A[i,1] = phi0(design[i,1], N)
    A[i,N+1] = phiN(design[i,1], N)
    for(j in 2 : (N)){
      A[i,j] = phii(design[i,1], j-1, N)
    }
  }
  
  Gamma <- matrix(data = 0, nrow = N+1, ncol = N+1)
  for(i in 1 : (N+1)){
    for(j in 1 : (N+1)){
      Gamma[i, j] = k(u[i], u[j], sig , theta)
    }
  }
  Gamma <- Gamma + nugget * diag(N+1)
  
  
  
  invGamma1 <- chol(Gamma)
  invGamma <- chol2inv(invGamma1)
  
  Amat2 <- diag((N+1))
  Amat1 <- rbind(A,Amat2)
  Amat <- rbind(Amat1,-Amat2)
  
  zetoil <- solve.QP(invGamma,dvec=rep(0,(N+1)),Amat=t(Amat),bvec=c(response,rep(lower,(N+1)),rep(-upper,(N+1))),meq=n)$solution
  
  return(structure(
    list(zetoil=zetoil,phi0=phi0,phiN=phiN,phii=phii,Amat=Amat,Gamma=Gamma, A=A, D=D,
         call=list(design=design,response=response,basis.size=basis.size,covtype=covtype,
                   coef.cov=coef.cov,coef.var=coef.var, nugget=nugget, lower=lower,
                   upper=upper)),class="kmBounded"))
}




#' @title Apply basis functions to given data in design space
#' @param model km* model
#' @param newdata data in design space
#' @method Phi1D kmBounded
Phi1D.kmBounded <- function(model, newdata){
  N <- model$call$basis.size  
  x <- newdata
  
  v <- matrix(0, nrow = length(x), ncol = N + 1)
  v[,1] <- apply(t(x), 2, model$phi0, N = N)
  v[,N+1] <- apply(t(x), 2, model$phiN, N = N)
  for(j in 2 : (N)){
    v[,j] = model$phii(x, j-1, N)
  }
  
  return(v)
}