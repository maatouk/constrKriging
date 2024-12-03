coef.cov_LOO <- function(object) {

  fctGamma <- object$fctGamma
  input <- object$call$design
  output <- object$call$response
  n <- length(output)
  
  fctinvGamma <- function(theta){
    invGamma1 <- chol(fctGamma(theta))
    invGamma <- chol2inv(invGamma1)
    return(invGamma)
  }
  
  Ai <- function(i) {
    object$A[-i, ]
  }
  
  Amat1 <- diag(nrow(object$Gamma))
  Amat1[1, 1] <- 0
  Amati <- function(i) {
    rbind(Ai(i), Amat1)
  }
  
  fctzetoil <- function(theta, i) {
    solve.QP(Dmat = fctinvGamma(theta), dvec = rep(0, nrow(object$Gamma)), Amat = t(Amati(i)), 
             bvec = c(output[-i], rep(0, nrow(object$Gamma))), meq = n-1)$solution
  }
  
  fctoptim <- function(theta) {
    sum <- 0
    for (i in 1 : length(input)) {
      sum = sum + (output[i] - (object$A%*%fctzetoil(theta, i))[i])^2
    }
    return(sum/length(input))
  }
  
  fctoptim.try <-  function(theta) {
    f <- NaN
    try(f <- fctoptim(theta))
    if(is.nan(f)) cat("(!) ",theta," -> ",f)
    return(f)
  }
  
  range_design <- range(input)
  lower <- min(dist(input)) * 1E-5
  upper <- (range_design[2]-range_design[1]) * 10
  
  #plot(Vectorize(fctoptim.try), xlim=c(lower,upper))
  
  ## 1-D optim, slow to converge, and often leads to high theta
  #tol <- (upper-lower)/1000
  #min <- optimize(fctoptim.try,lower=lower,upper = upper,tol = tol)$minimum
  #min <- optim(par = object$call$coef.cov,fctoptim.try,method = "L-BFGS-B",lower=lower,upper =upper,abstol=tol)$par
  
  ## search for a not(too-hight theta, which gives a reasonable LOO. Faster and more stable.
  thetas <- seq(from = lower, to = upper, length = 20)
  mins <- Vectorize(fctoptim.try)(thetas)
  if (all(is.nan(mins))) 
    stop("Error: Could not evaluate LOO at any theta !")
  #plot(thetas, mins)
  lower_bound_mins <- min(mins)+(max(mins)-min(mins))*.01
  lower_theta <- min(thetas[which(mins<lower_bound_mins)])
  
  return(lower_theta)
}
