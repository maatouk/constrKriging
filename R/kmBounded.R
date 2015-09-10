
#' @title Kriging model with bounds constraints
#' @param design 1-column matrix of the design of experiments
#' @param response a vector containing the output values given by the real function at the design points
#' @param basis.size a value represents the number of the basis functions (descritization of 1D input set)
#' @param covtype an optimal character string specifying the covariance function to be used ("gauss" and "matern3_2" choice)
#' @param coef.cov a value corresponding to the length theta hyper-parameters of covariance function
#' @param coef.var a value specifying the variance parameter
#' @param nugget an optimal value used as nugget effect to solve the numerical inverse matrix problem
#' @param lower lower bound constraint
#' @param upper upper bound constraint

#' @examples 
#' model = kmBounded1D(design=c(0.1, 0.3, 0.5, 0.9), response=c(7, -8, 9, 15), lower=-10, upper = 18)
kmBounded1D <- function(design, response, 
                        basis.size = dim(design)[1]+2+10, 
                        covtype = "gauss",
                        coef.cov = "LOO", #0.5*(max(design)-min(design)),
                        coef.var = var(response),
                        lower = min(response)-(max(response)-min(response))*.1,
                        upper = max(response)+(max(response)-min(response))*.1,
                        nugget = 1e-7*sd(response)) {
  
  stop("This model is under construction.")
  }
