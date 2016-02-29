#' @title Build a spline constrained function 
#' @param object km* model 
#' constrSpline(object=kmMonotonic1D(design=c(0.1, 0.2, 0.9), response=c(1, 5, 9)))(seq(f=0,t=1,l=100))
#' constrSpline(object=kmConvex1D(design=c(0.1, 0.5, 0.9), response=c(10, 5, 9)))(seq(f=0,t=1,l=100))
#' @examples 
#' design = c(0.05, 0.18, 0.3, 0.5, 0.7, 0.8)
#' response = c(-5, 2, 4, 5, 9, 15)
#' object = kmMonotonic1D(design, response)
#' f = constrSpline(object)
#' curve(f)
#' points(design, response, pch=19)
constrSpline <- function(object){
  f = function(x){
    Phi1D(model=object, newdata=x)%*%matrix(object$zetoil, ncol=1)
  }
  return(f)
}

#' @title Generic Phi1D.km* function
#' @param model km* model
#' @param newdata data in design space
Phi1D <- function(model, newdata) UseMethod("Phi1D")

