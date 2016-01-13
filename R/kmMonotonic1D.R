
#' @title Kriging model with monotonicity 1D constraint
#' @param design 1-column matrix of the design of experiments
#' @param response a vector containing the output values given by the real function at the design points
#' @param basis.size a value represents the number of the basis functions (descritization of 1D input set)
#' @param covtype an optimal character string specifying the covariance function to be used ("gauss" and "matern3_2" choice)
#' @param basis.type an optimal character string specifying the regularization of the basis functions to be used (which implies the class of the simulation paths)
#' @param coef.cov a value corresponding to the length theta hyper-parameters of covariance function
#' @param coef.var a value specifying the variance parameter
#' @param nugget an optimal value used as nugget effect to solve the numerical inverse matrix problem
#' @import quadrpog

#' @examples 
#' kmMonotonic1D(design=c(0.1, 0.5, 0.9), response=c(1, 5, 9))
#' kmMonotonic1D(design=c(0.1, 0.5, 0.9), response=c(1, 5, 9), basis.type="C0")
#' kmMonotonic1D(design=sort(runif(5)), response=cumsum(runif(5)), nugget=1e-4, covtype="gauss")

# ## Golchi Example
#' f <- function(x){
#' log(20*x+1)
#' }
#' design <- c(0, 0.1, 0.2, 0.3, 0.4, 0.9, 1)
#' response <- f(design)
#' meany <- mean(response)
#' f <- function(x){
#' log(20*x+1)-meany
#' }
#' design <- c(0, 0.1, 0.2, 0.3, 0.4, 0.9, 1)
#' response <- f(design)
#' model = kmMonotonic1D(design, response, coef.var=1, coef.cov=0.3,basis.size=50)



kmMonotonic1D <- function(design, response, 
                          basis.size = dim(design)[1]+2+10, 
                          covtype = "matern5_2",
                          basis.type = "C1", 
                          coef.cov = 0.5*(max(design)-min(design)),
                          coef.var = var(response),
                          nugget = 1e-7*sd(response)) {
  
  if (!is.matrix(design)) design=matrix(design,ncol=1)
  
  #   if (coef.cov=="LOO") {
  #     object=kmMonotonic1D(design, response, 
  #                          basis.size, 
  #                          covtype,
  #                          basis.type, 
  #                          coef.cov = 0.5*(max(design)-min(design)),
  #                          coef.var,
  #                          nugget)
  #     
  #     theta=coef.cov_LOO(object)
  #     
  #     model = NULL
  #     while(is.null(model)) # auto raise nugget if needed
  #       try(model <- kmMonotonic1D(design, response, 
  #                                  basis.size, 
  #                                  covtype ,
  #                                  basis.type,
  #                                  coef.cov = theta,
  #                                  coef.var,
  #                                  nugget=nugget*10))
  #     
  #     return(model)
  #   }
  
  
  n <- nrow(design) # number of design points
  N <- basis.size   # discretization size
  p <- (N + 2) - n  # degree of freedom
  u <- seq(0, 1, by = 1/N)  # discretization vector
  
  sig <- sqrt(coef.var)
  theta <- coef.cov
  
  
  if(covtype=='gauss'){
    # Gaussian covariance kernel
    k <- function(x, xp, sig, theta){
      (sig^2)*exp(-(x-xp)^2/(2*theta^2))
    }
    # derivative of the Gaussian kernel with respect to the first variable
    kp1 <- function(x, xp, sig, theta){
      -(x-xp)/(theta^2)*k(x,xp, sig, theta)
    }
    # derivative with respect to the second variable
    kp2 <- function(x, xp, sig, theta){
      -kp1(x,xp, sig, theta)
    }
    # second derivative with respect to the 1er and second variables
    kpp <- function(x, xp, sig, theta){
      (1/(theta^2))*k(x, xp, sig, theta)*(1-(x-xp)^2/(theta^2))
    }
    
  }
  else if (covtype=='matern3_2'){
    # Matern 3/2 covariance kernel
    k <- function(x, xp, sig, theta){
      (sig^2)*(1+(sqrt(3)*abs(x-xp)/theta))*exp(-sqrt(3)*abs(x-xp)/theta)
    }
    # derivative with respect to the first variable
    kp1 <- function(x, xp, sig, theta){
      sig^2*(sqrt(3)/theta*signp(x-xp)*exp(-sqrt(3)*(abs(x-xp)/theta))*(-sqrt(3)/theta*abs(x-xp)))
    }
    # derivative with respect to the second variable
    kp2 <- function(x, xp, sig, theta){
      -kp1(x,xp, sig, theta)
    }
    # second derivative with respect to the first and second variables
    kpp <- function(x, xp, sig, theta){
      sig^2*((3/theta^2)*exp(-sqrt(3)/theta*abs(x-xp))*(1-sqrt(3)/theta*abs(x-xp)))
    }
    
  }
  else if(covtype=='matern5_2'){    
    # Matern 5/2 covariance kernel
    k <- function(x, xp, sig, theta){
      sig^2*(1+sqrt(5)*(abs(x-xp))/theta+(5*(x-xp)^2)/(3*theta^2))*exp(-sqrt(5)*(abs(x-xp))/theta)
    }
    # derivative with respect to the first variable
    kp1 <- function(x, xp, sig, theta){
      sig^2*(sqrt(5)/(theta)*signp(x-xp)+10*(x-xp)/(3*theta^2))*exp(-sqrt(5)*(abs(x-xp))/theta)-
        sqrt(5)/(theta)*signp(x-xp)*k(x, xp, sig, theta)
    }
    # derivative with respect to the second variable
    kp2 <- function(x, xp, sig, theta){
      -kp1(x, xp, sig, theta)
    }
    # second derivative with respect to the first and second variables
    kpp <- function(x, xp, sig, theta){
      sig^2*(-10/(3*theta^2)*exp(-sqrt(5)*(abs(x-xp))/theta)+
               (sqrt(5)/theta)*signp(x-xp)*exp(-sqrt(5)*(abs(x-xp))/theta)*(sqrt(5)/theta*signp(x-xp)+
                                                                              10*(x-xp)/(3*theta^2)))-(sqrt(5)/theta)*signp(x-xp)*kp2(x, xp, sig, theta)
    }
  }
  
  else stop('covtype', covtype, 'not supported')
  
  
  
  if(basis.type=="C0"){
    ## basis functions (hat functions)      
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
      phi((x-u[i+1])*N)
    }
    
    A <- matrix(data = NA, ncol = N+1, nrow = n)
    A[, 1] = phi0(design[, 1], N)
    A[, N+1] = phiN(design[, 1], N)
    for(j in 2 : N){
      A[, j] = phii(design[, 1], j-1, N)
    }
    
    fctGamma=function(.theta){
      Gamma <- matrix(data = NA, nrow = N+1, ncol = N+1)
      for(i in 1 : (N+1)){
        for(j in 1 : (N+1)){
          Gamma[i, j] = k(u[i], u[j], sig , .theta)
        }
      }
      Gamma <- Gamma + nugget * diag(N+1)
      return(Gamma)
    }
    Gamma=fctGamma(theta)
    
    
  }else if (basis.type == "C1"){
    ## basis functions (hat functions)
    h <- function(x){
      ifelse(x >= -1 & x <= 1, 1-abs(x), 0)
    }
    h0 <- function(x, N){
      h(x*N)
    }
    hN <- function(x, N){
      h((x-u[N+1])*N)
    }
    hi <- function(x, i, N){
      h((x - u[i+1])*N)
    }
    ## primitive of the hat functions
    phii <- function(x, i, N){
      delta <- 1/N
      ifelse(x <= u[i], 0,  
             ifelse(x >= u[i] & x <= u[i+1], ((x-u[i])*hi(x, i, N))/2,
                    ifelse(x >= u[i+1] & x <= u[i+2], delta - ((u[i+2]-x)*hi(x, i, N))/2,
                           delta
                    )
             )
      )
    }
    phi0 <- function(x, N){
      delta <- 1/N
      ifelse(x <= -delta, -delta/2,  
             ifelse(x >= -delta & x <= 0, -delta/2+((x+delta)*h0(x, N))/2,
                    ifelse(x >= 0 & x <= delta, delta/2 - ((delta-x)*h0(x, N))/2, delta/2
                    )
             )
      )
    }
    phiN <- function(x, N){
      delta <- 1/N
      ifelse(x >= u[N] & x <= u[N+1], (x - u[N])*hN(x, N)/2, 
             ifelse(x<= u[N],0, ifelse(x>=u[N+1] & x<= 1+delta, 
                                       delta-(1+delta-x)*hN(x, N)/2, delta)))
    }
    
    fctGamma=function(.theta){
      Gamma <- matrix(data = NA, nrow = N+2, ncol = N+2)
      Gamma[1,1] <- k(0,0, sig, .theta)
      for(j in 2 : (N+2)){
        Gamma[1,j] <- kp2(0, u[j-1], sig, .theta)
      }
      for(i in 2 : (N+2)){
        Gamma[i,1] <- kp1(u[i-1], 0, sig, .theta)
      }
      for(i in 2 : (N+2)){
        for(j in 2 : (N+2)){
          Gamma[i, j] = kpp(u[i-1], u[j-1], sig, .theta)
        }
      }
      Gamma <- Gamma + nugget * diag(N+2) 
      return(Gamma)
    }
    Gamma=fctGamma(theta)
    
    A <- matrix(data = NA, ncol = N+2, nrow = n)
    A[, 1] = 1
    A[, 2] = phi0(design[, 1], N)
    A[, N+2] = phiN(design[, 1], N)
    for(j in 3 : (N+1)){
      A[, j] = phii(design[, 1], j-2, N)
    }
  }
  
  else stop ("basis.type",basis.type, "not supported")
  
  invGamma <- chol2inv(chol(Gamma))
  
  if (basis.type=="C0"){
    CN <- diag(N)
    for(i in 1 : N){
      for(j in 1 : N){
        if(j==i+1){
          CN[i, j] = -1
        }
      }
    }
    v <- rep(0, N)
    v[N] <- -1
    
    BN <- cbind(CN, v)
    D <- -BN
    Amat <- rbind(A, D)
  }
  else {
    Amat1 <- diag(N+2)
    Amat1[1, 1] <- 0
    Amat <- rbind(A, Amat1)
  }
  
  
  if (basis.type=="C0")
    zetoil <- solve.QP(invGamma, dvec=rep(0, N+1), Amat=t(Amat), bvec=c(response, rep(0, N)), meq=n)$solution
  else
    zetoil <- solve.QP(invGamma, dvec=rep(0, N+2), Amat=t(Amat), bvec=c(response, rep(0, N+2)), meq=n)$solution
  
  return(structure(
    list(zetoil=zetoil, phi0=phi0, phiN=phiN, phii=phii, Amat=Amat, Gamma=Gamma, A=A, D=D, fctGamma=fctGamma, invGamma=invGamma,
         call=list(design=design, response=response, basis.size=basis.size, covtype=covtype, basis.type=basis.type,
                   coef.cov=coef.cov, coef.var=coef.var, nugget=nugget)
    ), class="kmMonotonic1D"))
}




#' @title Apply basis functions to given data in design space
#' @param model km* model
#' @param newdata data in design space
#' @method Phi1D kmMonotonic1D
Phi1D.kmMonotonic1D <- function(model, newdata){
  N <- model$call$basis.size  
  x <- newdata
  if(model$call$basis.type == 'C0'){
    v <- matrix(NA, nrow = length(x), ncol = N + 1)
    v[, 1] <- model$phi0(x, N = N)
    v[, N+1] <- model$phiN(x, N = N)
    for(j in 2 : N){
      v[, j] = model$phii(x, j-1, N)
    }
    
    
  }else {
    
    v <- matrix(NA, nrow = length(x), ncol = N + 2)
    v[, 1] <- 1
    v[, 2] <- model$phi0(x, N = N)
    v[, N+2] <- model$phiN(x, N = N)
    for(j in 3 : (N+1)){
      v[, j] = model$phii(x, j-2, N)
    }
  }
  return(v)
}

#' @test signp(rnorm(10))
#' @test signp(rpois(100,1))
signp <- function(t){
  s=sign(t)
  s[which(s==0)]=1
  s
}