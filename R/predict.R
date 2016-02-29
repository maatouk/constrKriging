#' @title Predict response from a km* model
#' @param object km* model 
#' @param newdata a vector which represents the points where to performs predictions
#' @param nsim the number of response vectors to simulate

#' @examples 
#' predict(object=kmMonotonic1D(design=c(0.1, 0.5, 0.9), response=c(1, 5, 9)), newdata=seq(f=0, t=1, l=100), nsim=10)
#' predict(object=kmConvex1D(design=c(0.1, 0.5, 0.9), response=c(5, 1, 9)), newdata=seq(f=0, t=1, l=100), nsim=10)

#' design = c(0.1, 0.5, 0.9)
#' response = c(1, 8, 9)
#' model = kmMonotonic1D(design, response)
#' graphics::plot(x=seq(0,1,,100),y=predict(object=model, newdata=seq(0,1,,100), nsim=100)[,'Median'], type='l')
#' points(design, response, pch=19)

predictCK <- function(object, newdata, nsim = 100){
  y <- simulate_process(object, nsim, seed=1, newdata)
  return(t(sapply(apply(y, ecdf, MARGIN = 1), summary)))
}

#' @title Generic simulate_process method
#' @param object km* model
#' @param nsim number of simulated processes
#' @param seed random seed
#' @param newdata data in design space where to evaluate simlations
simulate_process <- function(object, nsim, seed, newdata) UseMethod("simulate_process")

predict.kmMonotonic1D <- predict.kmConvex1D <- predict.kmBounded1D <- predictCK