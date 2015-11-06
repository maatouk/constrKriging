#' @title Simulate responses vectors from a kmMonotonic1D model
#' @param object kmMonotonic1D model
#' @param newdata a vector which represents the points where to performs predictions
#' @param nsim the number of response vectors to simulate
#' @param seed optional random seed
#' @import MASS

#' simulate(object=kmMonotonic1D(design=c(0.1, 0.5, 0.9), response=c(1, 5, 9)),newdata=seq(f=0,t=1,l=100),nsim=50)
#' @examples 
#' design = c(0.1, 0.5, 0.9)
#' response = c(1, 5, 5.5)
#' model <- kmMonotonic1D(design, response, covtype="matern5_2", coef.var=1, basis.size=50)
#' x=seq(0,1,,100)
#' graphics::matplot(x,y=simulate_process(object=model, newdata=x, nsim=40),col='gray',type='l',lty=1)
#' lines(x,constrSpline(object=model)(x), lty=1,col='black')
#' points(design, response, pch=19)


#' ## Golchi Example
#' f <- function(x){
#' log(20*x+1)
#' }
#' design <- c(0, 0.1, 0.2, 0.3, 0.4, 0.9, 1)
#' response <- f(design)
#' meany <- mean(response)
#' f <- function(x){
#'  log(20*x+1)-meany
#' }
#' design <- c(0, 0.1, 0.2, 0.3, 0.4, 0.9, 1)
#' response <- f(design)
#' model = kmMonotonic1D(design, response, covtype="matern5_2", coef.var=50, coef.cov=2.5, basis.size=50)
#' graphics::matplot(x=seq(0,1,,100),y=simulate_process(object=model, newdata=seq(0,1,,100), nsim=40),col='gray',type='l',lty=1)
#' lines(x,constrSpline(model)(x),lty=2,col='black',lwd=1)
#' lines(x,f(x))
#' points(design,response,pch=19)

simulate_process.kmMonotonic1D <- function(object, nsim, seed=NULL, newdata){
  if (!is.null(seed)) set.seed(seed)
  
  N <- object$call$basis.size  
  zetoil <- object$zetoil
  A <- object$A
  Gamma <- object$Gamma
  p <- ncol(A)-nrow(object$call$design)
  response <- object$call$response
  D <- object$D
  
  
  B <- diag(ncol(A)) - t(A) %*% chol2inv(chol(A %*% t(A))) %*% A     
  epsilontilde <- eigen(t(B) %*% chol2inv(chol(Gamma)) %*% B)$vectors
  epsilon <- B %*% epsilontilde[, 1 : p]
  c <- eigen(t(B) %*% chol2inv(chol(Gamma)) %*% B)$values[1:p]
  d <- 1/c[1 : p]                     
  
  
  zcentre <- Gamma %*% t(A) %*% chol2inv(chol(A %*% Gamma %*% t(A))) %*% response
  
  
  setoil <- t(epsilon) %*% (zetoil - zcentre)
  
if (object$call$basis.type == 'C0'){
    Xi <- matrix(-1, ncol = nsim, nrow = ((ncol(A))))
    for (j in 1 : nsim){
      Xi_current <- Xi[, j]
      unif <- 1
      t <- 0
      while(unif > t ){
        Xi_current <- Xi[, j]
        while (min(D%*%Xi_current) < 0){
          #    s <- rnorm(1, setoil, d) ## Dans le Cas o? p = 1 i.e (N = n). 
#           s <- matrix(mvrnorm(1, as.vector(setoil), diag(d)), ncol = 1)
          s <- setoil + sqrt(d)*matrix(rnorm(p, 0, 1), ncol = 1)
          Xi_current <- as.vector(zcentre) + (epsilon %*% s)
        }
        t <- as.numeric(exp(sum((setoil -s) * setoil * c)))  
        unif <- runif(1, 0, 1)
      }
      Xi[, j] <- Xi_current
    } 
  } else {
    Xi <- matrix(-1, ncol = nsim, nrow = ncol(A))
    for (j in 1 : nsim){
      Xi_current <- Xi[, j]
      unif <- 1
      t <- 0
      while(unif > t ){
        while (min(Xi_current[-1]) < 0){
          #s <- rnorm(1, 0, d) ## Dans le Cas o? p = 1 i.e (N = n). 
#           s <- matrix(mvrnorm(1, as.vector(setoil), diag(d)), ncol = 1)
          s <- setoil + sqrt(d)*matrix(rnorm(p, 0, 1), ncol = 1)
          Xi_current <- as.vector(zcentre) + (epsilon %*% s)
        }
        t <- as.numeric(exp(sum(setoil * setoil / d) - sum(s * setoil / d)))
        unif <- runif(1, 0, 1)
      }
      Xi[, j] <- Xi_current
    }
  }
  
  
  return(Phi1D.kmMonotonic1D(object, newdata)%*%Xi)
}