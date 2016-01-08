#' @title Simulate responses vectors from a kmBounded model
#' @param object kmBounded model
#' @param newdata a vector which represents the points where to performs predictions
#' @param nsim the number of response vectors to simulate
#' @param seed optional random seed
#' @import MASS
#' @import solve.QP

#' @examples 
#' design=c(0, 0.1, 0.2, 0.42, 0.5, 0.9)
#' response = c(10, 7, -8, -5, 10, 15)
#' model = kmBounded1D(design, response, lower = -20, upper = 20, coef.cov=0.2, coef.var=100, basis.size = 50)
#' x = seq(0,1,,101)
#' graphics::matplot(x,y=simulate_process(object=model, newdata=x, nsim=100),
#' type='l', col='gray', lty = 1, ylab = "response", ylim=c(model$call$lower,model$call$upper))
#' lines(x,constrSpline(object=model)(x), lty=1,col='black')
#' points(design, response, pch=19)
#' abline(h=model$call$lower, lty=2)
#' abline(h=model$call$upper, lty=2)

simulate_process.kmBounded1D <- function(object, nsim, seed=NULL, newdata){
  if (!is.null(seed)) set.seed(seed)
  
  N <- object$call$basis.size  
  zetoil <- object$zetoil
  A <- object$A
  Gamma <- object$Gamma
  invGamma <- object$invGamma
  p <- ncol(A)-nrow(object$call$design)
  response <- object$call$response
  D <- object$D
  lower <- object$call$lower
  upper <- object$call$upper
  
  B <- diag(ncol(A)) - t(A) %*% chol2inv(chol(A %*% t(A))) %*% A    
  epsilontilde <- eigen(t(B) %*% invGamma %*% B)$vectors
  epsilon <- B %*% epsilontilde[, 1 : p]
  c <- eigen(t(B) %*% invGamma %*% B)$values[1 : p]
  d <- 1 / c[1 : p]       
  zcentre <- Gamma %*% t(A) %*% chol2inv(chol(A %*% Gamma %*% t(A))) %*% response 
  setoil <- t(epsilon) %*% (zetoil - zcentre)
  
  Xi <- matrix(-10+lower, ncol = nsim, nrow = N+1)
  for (j in 1 : nsim){
    Xi_current <- Xi[, j]
    unif <- 1
    t <- 0  
    while(unif > t){
      Xi_current <- Xi[, j]
      while ((min(Xi_current) < lower) || (max(Xi_current) > upper)){
        s <- setoil + sqrt(d)*matrix(rnorm(p, 0, 1), ncol = 1)
        Xi_current <- as.vector(zcentre) + (epsilon %*% s)
      }
      t <- as.numeric(exp(sum((setoil - s) * setoil * c)))      
      unif <- runif(1, 0, 1)
    }
    Xi[, j] <- Xi_current
  }
  
  
  return(Phi1D.kmBounded1D(object, newdata)%*%Xi)
}