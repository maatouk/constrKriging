coef.cov_LOO =function(object) {
  
  fctGamma = object$fctGamma
  
  input=object$call$design
  output=object$call$response
  
  n=length(output)
  
  fctinvGamma <- function(theta){
    invGamma1 = chol(fctGamma(theta))
    invGamma = chol2inv(invGamma1)
    return(invGamma)
  }
  
  
  Ai <- function(i){
    object$A[-i,]
  }
  
  
  Amat1 <- diag(nrow(object$Gamma))
  Amat1[1,1] <- 0
  #Amat <- rbind(object$A,Amat1)
  Amati <- function(i){
    rbind(Ai(i),Amat1)
  }
  
  fctzetoil <- function(theta,i){
    solve.QP(fctinvGamma(theta),dvec=rep(0,nrow(object$Gamma)),Amat=t(Amati(i)),bvec=c(output[-i],rep(0,nrow(object$Gamma))),meq=n-1)$solution
  }
  
  
  fctoptim <- function(theta){
    sum = 0
    for(i in 1 : length(input)){
      sum = sum + (output[i]-(object$A%*%fctzetoil(theta,i))[i])^2
    }
    return(sum)
  }
  
  range_design=range(input) 
#plot(Vectorize(fctoptim), xlim=c(min(dist(input))*1E-5,(range_design[2]-range_design[1])*2))
  optim(par = object$call$coef.cov,fctoptim,method = "L-BFGS-B",lower=min(dist(input))*1E-5,upper = (range_design[2]-range_design[1])*2)
}