#' @title Simulate responses vectors from a kmBounded model
#' @param object kmBounded model
#' @param newdata a vector which represents the points where to performs predictions
#' @param nsim the number of response vectors to simulate
#' @param seed optional random seed
#' @import MASS
#' 
#' simulate_process(object=kmBounded(design=c(0.1, 0.3, 0.5, 0.9), response=c(7, -8, 9, 15)),newdata=seq(f=0,t=1,l=100),nsim=10)
#' @examples 
#' design=c(0.1, 0.3, 0.5, 0.9)
#' response = c(7, -8, 9, 15)
#' model = kmBounded(design, response, lower = -10, upper = 18)
#' x=seq(0,1,,100)
#' graphics::plot(x=x,y=simulate_process(object=model, newdata=x, nsim=1),type='l')

simulate_process.kmBounded <- function(object, nsim, seed=NULL, newdata){
  if (!is.null(seed)) set.seed(seed)
  
  N <- object$call$basis.size  
  zetoil <- object$zetoil
  A <- object$A
  Gamma <- object$Gamma
  p <- ncol(A)-nrow(object$call$design)
  response <- object$call$response
  D <- object$D
  lower <- object$call$lower
  upper <- object$call$upper
  
  B <- diag(ncol(A)) - t(A) %*% solve(A %*% t(A)) %*% A     ## Proj sur F_0
  epsilontilde <- eigen(t(B) %*% Gamma %*% B)$vectors
  epsilon <- B %*% epsilontilde[, 1 : p]
  c <- eigen(t(B) %*% Gamma %*% B)$values
  d <- c[1 : p]                     # le plus grand p valeurs propres
  
  
  zcentre <- Gamma %*% t(A) %*% solve(A %*% Gamma %*% t(A)) %*% response
  
  
  setoil <- t(epsilon) %*% (zetoil - zcentre)
  
  Xi <- matrix(-10+lower, ncol = nsim, nrow = N+1)
  for (j in 1 : nsim){
    Xi_current <- Xi[,j]
    while ((min(Xi_current) < lower) || (max(Xi_current) > upper)){
      #s <- rnorm(1, 0, d) ## Dans le Cas o? p = 1 i.e (N = n). 
      s <- matrix(mvrnorm(1, as.vector(setoil), diag(d)), ncol = 1)
      t <- as.numeric(exp(sum(setoil * setoil / d) - sum(s * setoil / d)))
      unif <- runif(1, 0, 1)
      while(unif > t ){
        #s <- rnorm(1, setoil, d) # Pour N = n
        s <- matrix(mvrnorm(1, as.vector(setoil), diag(d)), ncol = 1)
        t <- as.numeric(exp(sum(setoil * setoil / d) - sum(s * setoil / d)))
        unif <- runif(1, 0, 1)
      }
      if((min(zcentre >= lower)) & (max(zcentre) <= upper)){
        s <- matrix(mvrnorm(1, rep(0, p), diag(d)))
        Xi_current <- as.vector(zcentre) + (epsilon %*% s)
      }
      else {
        Xi_current <- as.vector(zetoil) + (epsilon %*% s)
      }
    }
    Xi[, j] <- Xi_current
  }
  
  
  return(Phi1D.kmBounded(object, newdata)%*%Xi)
}