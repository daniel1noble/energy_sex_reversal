# Functions

#' @title pmcmc
#' @description Calculates the, p-value or pMCMC value for a posterior distribution
#' @param x The vector for the posterior distribution. Note that this will test the null hypothesis that the parameter of interest is significantly different from 0. 

pmcmc <- function(x){
  2*(1 - max(table(x<0) / nrow(x)))
}