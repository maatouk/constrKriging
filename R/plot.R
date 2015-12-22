#' @title Plot the km* object
#' @param x The coordinates to plot in the object design space (seq(0,1,,100) by default)
#' @param object The km* model to plot
#' @param spline Whether to plot or not the spline of the model (by default TRUE)
#' @param nsim Number of simulations of the model realizations to build the sample in order to get the median, mean, minmax and quantiles statistics (by default 100)
#' @param median Whether to plot or not the median of the model sample (by default TRUE)
#' @param mean Whether to plot or not the mean of the model sample (by default FALSE)
#' @param minmax Whether to plot or not the min/max of the model sample (by default FALSE)
#' @param quantiles Whether to plot or not the 1st and 3rd quantiles of the model sample (by default TRUE)
#' @param col Color to use for drawing
#' @param add Should we add this plot to an existing one ?
#' @export

#' @examples 
#' ## Convex Example
#' design = c(0.1, 0.5, 0.9)
#' response = c(10, 5, 9)
#' model = kmConvex1D(design, response)
#' plot(object=model, spline=TRUE, quantiles=TRUE, minmax=FALSE, col='gray',nsim=10)
#' points(design,response,pch=19)

#' ## Monotone Example
#' design = c(0.1, 0.5, 0.9)
#' response = c(1, 5, 9)
#' model = kmMonotonic1D(design, response, covtype="matern5_2", coef.cov=1, coef.var=100)
#' plot(object=model, spline=TRUE, median=FALSE, quantiles=TRUE, minmax=FALSE, col='gray',nsim=100)
#' points(design,response,pch=19)

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
#' model = kmMonotonic1D(design, response, covtype="matern5_2", coef.var=100, coef.cov=2, basis.size=50)
#' plot(object=model, median=FALSE, spline=TRUE, quantiles=TRUE, minmax=FALSE, col='gray',nsim=100)
#' points(design,response,pch=19)

#' ## Bounded Example
#' design <- c(0.1, 0.3, 0.5, 0.9)
#' response <- c(7, -8, 9, 15)
#' model = kmBounded1D(design, response, lower=-10, upper = 18, coef.cov=0.2, basis.size=50)
#' plot(object=model, median=FALSE, spline=TRUE, quantiles=TRUE, minmax=FALSE, col='gray',nsim=100)
#' points(design,response,pch=19)

plotCK <- function(x=seq(f=min(object$call$design),t=max(object$call$design),l=100), object, spline=FALSE, nsim=100, median=FALSE, mean=FALSE, minmax=FALSE, quantiles=TRUE, col='black',add=F){
  if (!isTRUE(add))
    graphics::plot(x=range(object$call$design), y=range(object$call$response), type='n', xlab='design', ylab='response')
  
  if (isTRUE(spline))
    lines(x,constrSpline(object)(x),lty=2,col=col)
  
  if (isTRUE(median) || isTRUE(mean) || isTRUE(quantiles) || isTRUE(minmax)) {
    pred <- predict(object, newdata=x, nsim )
    if (isTRUE(median))
      lines(x,pred[,'Median'],lty=2,col=col)
    
    if (isTRUE(mean))
      lines(x,pred[,'Mean'],lty=3,col=col)
    
    if (isTRUE(quantiles)){
      col.rgb=col2rgb(col)
      polygon(x=c(x,rev(x)), y=c(pred[,'1st Qu.'], rev(pred[,'3rd Qu.'])) , border = F,  col = rgb(col.rgb[1]/255,col.rgb[2]/255,col.rgb[3]/255,alpha = .4))
    }
    if (isTRUE(minmax)){
      col.rgb=col2rgb(col)
      polygon(x=c(x,rev(x)), y=c(pred[,'Min.'], rev(pred[,'Max.'])) ,  border = F,    col = rgb(col.rgb[1]/255,col.rgb[2]/255,col.rgb[3]/255,alpha = .2))
    }
  }
}

plot.kmMonotonic1D <- plot.kmConvex1D <- plot.kmBounded1D <- plotCK