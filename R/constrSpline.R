#' @title Build a spline constrained function 
#' @param object km* model 
#' constrSpline(object=kmMonotonic1D(design=c(0.1, 0.2, 0.9), response=c(1, 5, 9)))(seq(f=0,t=1,l=100))
#' constrSpline(object=kmConvex1D(design=c(0.1, 0.5, 0.9), response=c(1, 5, 9)))(seq(f=0,t=1,l=100))
#' @examples 
#' design = c(0.1, 0.5, 0.9)
#' response = c(1, 5, 9)
#' object = kmMonotonic1D(design, response, basis.type="C1")
#' f = constrSpline(object)
#' curve(f)
constrSpline <- function(object){
  f = function(x){
    Phi1D(model=object, newdata=x)%*%matrix(object$zetoil,ncol=1)
  }
  return(f)
}

#' @title Generic Phi1D.km* function
#' @param model km* model
#' @param newdata data in design space
Phi1D <- function(model, newdata) UseMethod("Phi1D")

