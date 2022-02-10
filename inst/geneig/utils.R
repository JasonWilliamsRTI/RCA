#' Simulate example data
#'
#' Simulate MVN example data from the \code{corr} data, which are a correlation
#' matrix
#' 
#' @param An 11x11 correlation matrix
#' 
#' @N Integer. The number of observations to simulation.
#' 
#' @seed Integer. A random seed for replication
#' 
#' @import MASS
#' 
#' @export
corrSim <- function(corr, N=1e6, seed=123)
{
  varNames <- c("information",
                "digitspan",
                "vocab",
                "arith",
                "comp",
                "simil",
                "piccomp",
                "picarr",
                "block",
                "objass",
                "digsym")
  
  set.seed(seed)
  
  X <- mvrnorm(N, rep(0, length(varNames)), corr)
  X <- data.frame(X)
  names(X) <- varNames
  
  X
}
