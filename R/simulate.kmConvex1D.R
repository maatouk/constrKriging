#' @title Simulate responses vectors from a kmConvex1D object
#' @param object kmConvex1D model
#' @param newdata a vector which represents the points where to performs predictions
#' @param nsim the number of response vectors to simulate
#' @param seed optional random seed
#' @import MASS

#' @examples 
#' design = c(0.1, 0.5, 0.9)
#' response = c(10, 5, 9)
#' model = kmConvex1D(design, response, coef.cov = 0.35)
#' x = seq(0, 1,, 100)
#' graphics::matplot(x, y=simulate_process(object=model, newdata=x, nsim=100), type='l', col='gray', lty=1, ylab='response')
#' lines(x,constrSpline(object=model)(x), lty=1, col='black')
#' points(design, response, pch=19)
#' legend(0.15, 12.5, c("convex GP sample paths", "posterior max"), 
#'        col = c('gray', 'black'), text.col = "black",
#'        lty = c(1, 1), pch=c(NA_integer_, NA_integer_), lwd = c(2, 2), text.font=1, box.lty=0, cex=1)


simulate_process.kmConvex1D <- function(object, nsim, seed=NULL, newdata){
  if (!is.null(seed)) set.seed(seed)
  
  N <- object$call$basis.size  
  zetoil <- object$zetoil
  A <- object$A
  Gamma <- object$Gamma
  invGamma <- object$invGamma
  p <- ncol(A) - nrow(object$call$design)
  response <- object$call$response
  D <- object$D
  
  
  B <- diag(ncol(A)) - t(A) %*% chol2inv(chol(A %*% t(A))) %*% A     
  epsilontilde <- eigen(t(B) %*% invGamma %*% B)$vectors
  epsilon <- B %*% epsilontilde[, 1 : p]
  c <- eigen(t(B) %*% invGamma %*% B)$values[1 : p]
  d <- 1/c[1 : p]                     
  zcentre <- Gamma %*% t(A) %*% chol2inv(chol(A %*% Gamma %*% t(A))) %*% response
  setoil <- t(epsilon) %*% (zetoil - zcentre)
  
  Xi <- matrix(-1, ncol = nsim, nrow = (N+3))
  for (j in 1 : nsim){
    Xi_current <- Xi[, j]
    unif <- 1
    t <- 0
    while(unif > t ){
      Xi_current <- Xi[, j]
      while (min(Xi_current[-(1:2)]) < 0){
        s <- setoil + sqrt(d)*matrix(rnorm(p, 0, 1), ncol = 1)
        Xi_current <- as.vector(zcentre) + (epsilon %*% s)
      }
      t <- as.numeric(exp(sum((setoil - s) * setoil * c)))  
      unif <- runif(1, 0, 1)
    }
    Xi[, j] <- Xi_current
  }
  
  
  return(Phi1D.kmConvex1D(object, newdata)%*%Xi)
}