#' @title Kriging model with boundedness constraints
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
#' model = kmBounded1D(design=c(0.1, 0.3, 0.5, 0.9), response=c(7, -8, 9, 15), lower=-10, upper = 18, coef.cov=1)
kmBounded1D <- function(design, response, 
                        basis.size = dim(design)[1]+2+10, 
                        covtype = "matern5_2",
                        coef.cov = 0.5*(max(design)-min(design)), # "LOO"
                        coef.var = var(response),
                        lower = min(response)-(max(response)-min(response))*.1,
                        upper = max(response)+(max(response)-min(response))*.1,
                        nugget = 1e-7*sd(response)) {
  
  if (!is.matrix(design)) design = matrix(design, ncol = 1)
  
  
  #   if (coef.cov == "LOO") {
  #     object = kmBounded1D(design, response, 
  #                                      basis.size = dim(design)[1]+2+10, 
  #                                      covtype,
  #                                      coef.cov=0.5*(max(design)-min(design)),
  #                                      coef.var,
  #                                      lower,
  #                                      upper,
  #                                      nugget)
  #       
  #       theta <- coef.cov_LOO(object)
  #     
  #     model <- NULL
  #     while(is.null(model)) # auto raise nugget if needed
  #       try(model <- kmBounded1D(design, response, 
  #                                    basis.size, 
  #                                    covtype,
  #                                    coef.cov = theta,
  #                                    coef.var,
  #                                    lower,
  #                                    upper,
  #                                    nugget))
  #     
  #     return(model)
  #   }
  
  n <- nrow(design) # number of design points
  N <- basis.size   # discretization size
  p <- (N + 1) - n  # degree of freedom
  u <- seq(from = 0, to = 1, by = 1/N) # discretization vector
  
  sig <- sqrt(coef.var)
  theta <- coef.cov
  
  
  if (covtype == 'gauss') {
    # Gaussian covariance kernel 
    k <- function(x, xp, sig, theta){
      sig^2*exp(-(x-xp)^2/(2*theta^2))
    }
  }else if (covtype == 'matern3_2') {
    # Matern 3/2 covariance kernel
    k <- function(x, xp, sig, theta) {
      sig^2*(1+(sqrt(3)*abs(x-xp)/theta))*exp(-sqrt(3)*abs(x-xp)/theta)
    }
  }
  else if (covtype == "matern5_2") {
    # Matern 5/2 covariance kernel
    k <- function(x, xp, sig, theta) {
      sig^2*(1+sqrt(5)*(abs(x-xp))/theta+(5*(x-xp)^2)/(3*theta^2))*exp(-sqrt(5)*(abs(x-xp))/theta)
    }
  }
  
  else stop('covtype', covtype, 'not supported')
  
  
  # Basis functions  
  phi <- function(x) {
    ifelse (x >= -1 & x <= 1, 1-abs(x), 0)
  }
  phi0 <- function(x, N) {
    phi(x*N)
  }
  phiN <- function(x, N) {
    phi((x-u[N+1])*N)
  }
  phii <- function(x, i, N) {
    phi((x-u[i+1])*N)
  }
  
  A <- matrix(data = NA, ncol = N+1, nrow = n)
  A[, 1] = phi0(design[, 1], N)
  A[, N+1] = phiN(design[, 1], N)
  for (j in 2 : N) {
    A[, j] <- phii(design[, 1], j-1, N)
  }
  
  fctGamma <- function(.theta) {
    Gamma <- matrix(data = NA, nrow = N+1, ncol = N+1)
    for (i in 1 : (N+1)) {
      for (j in 1 : (N+1)) {
        Gamma[i, j] = k(u[i], u[j], sig, .theta)
      }
    }
    Gamma <- Gamma + nugget * diag(N+1)
    return(Gamma)
  }
  Gamma <- fctGamma(theta)
  
  invGamma <- chol2inv(chol(Gamma))
  
  Amat2 <- diag(N+1)
  Amat1 <- rbind(A, Amat2)
  Amat <- rbind(Amat1, -Amat2)
  
  zetoil <- solve.QP(Dmat = invGamma, dvec = rep(0, N+1), Amat = t(Amat), 
                     bvec = c(response, rep(lower, N+1), rep(-upper, N+1)), meq = n)$solution
  
  return(structure(
    list(zetoil = zetoil, phi0 = phi0, phiN = phiN, phii = phii, Amat = Amat, Gamma = Gamma, A = A, D = D, 
         fctGamma = fctGamma, invGamma = invGamma,
         call = list(design = design, response = response, basis.size = basis.size, covtype = covtype,
                   coef.cov = coef.cov, coef.var = coef.var, nugget = nugget, lower = lower,
                   upper = upper)), class = "kmBounded1D"))
}




#' @title Apply basis functions to given data in design space
#' @param model km* model
#' @param newdata data in design space
#' @method Phi1D kmBounded
Phi1D.kmBounded1D <- function(model, newdata) {
  N <- model$call$basis.size  
  x <- newdata
  
  v <- matrix(data = NA, nrow = length(x), ncol = N+1)
  v[, 1] <- model$phi0(x, N = N)
  v[, N+1] <- model$phiN(x, N = N)
  for (j in 2 : N) {
    v[, j] <- model$phii(x, j-1, N = N)
  }
  
  return(v)
}
