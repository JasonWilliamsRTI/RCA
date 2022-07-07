#' Reliable Component Analysis of Cliff & Caruso 1998
#' 
#' R implementation of Reliable Component Analysis of Cliff & Caruso 1998
#' 
#'   
#' @param data data frame with variables in \code{varlist}.
#' 
#' @param varlist Character list with the variable names in \code{data} that
#' should be included in the analysis.
#' 
#' @param reliab Numeric vector of reliabilities in (0,1)
#'
#' @param retain Integer. The number of components to retain.
#' 
#' @param criterion Value in (0,1), sets the minimum reliability
#'   
#' @param method See \code{\link{cor}}.
#'  
#' @param use See \code{\link{cor}}.
#' 
#' @param file Character. File name with an'.txt' extention for saving output.
#' 
#' @import geigen
#' @export
#' 
#' @examples 
#' 
#' # Replicate the analyses of Cliff and Caruso (1998)
#' 
#' library(RCA)
#' 
#' # simulate data using the Cliff and Caruso (1998) correlation matrix
#' data(corr)
#' data <- corrSim(corr)
#' reliab <- c(.90, .89, .96, .86, .84, .82, .85, .76, .89, .68, .86)
#' 
#' # view the RelComp documentation
#' ?RelComp
#' 
#' # run RelComp on the simulated data
#' output <- RelComp(data, names(data), reliab)
#' 
#' # view the reliabilities, compare to the first row in Table 8 of 
#' # Cliff and Caruso (1998)
#' round(output$reliabilities[1:5], 3)
#' 

RelComp <- function(data, varlist, reliab, 
                    retain=5, criterion=.70, 
                    method = "pearson", 
                    use = "pairwise.complete.obs",
                    file = "RelCompOutput.txt")
{
  
  # estimate the correlation matrix from the data
  corr <- cor(data[, varlist], method = method, use = use)
  
  # create rstar matrix with corrs on off-diagonal and reliabilitys on diagonal
  diag  <- diag(diag(corr))
  rel2  <- diag(reliab)
  diff  <- diag - rel2
  rstar <- corr - diff
  
  # generalized eigenvalues - equal to component reliabilities
  # generalized eigenvectors - vectors of weights that are not
  # normed to have sum of squared elements = 1
  reliabilities   <- geigen(rstar, corr)
  weights_nonorm  <- reliabilities$vectors
  weights_nonorm  <- -1*weights_nonorm[,ncol(weights_nonorm):1] #TODO: jason sending output, why columns not all -1*?
  reliabilities   <- rev( reliabilities$values )
  
  nvar <- nrow(weights_nonorm)
  
  # normed weights
  weights_normed <- weights_nonorm
  for(i in 1:nvar)
  {
    weights_normed[,i] <- weights_normed[,i]/norm(t(as.matrix(weights_nonorm[,i])))
  }
  
  # subset weights matrices to retained compoents 
  # get loadings by multiplying (de normed) weights by correlation matrix
  
  # create vector of cutoff indicators
  critcomp <- reliabilities >= criterion
  critrel  <- sum(critcomp)
  
  # specify ncomps as retain and superseded if critrel > 0
  ncomps <- retain 
  if(critrel > 0) ncomps <- critrel
  
  # subset data 
  weights_nonorm_sub <- weights_nonorm[, 1:ncomps]
  weights_normed_sub <- weights_normed[, 1:ncomps]
  # QUESTION: shouldn't we select columns based on critcomp?
  #weights_nonorm_sub <- weights_nonorm[, critcomp]
  #weights_normed_sub <- weights_normed[, critcomp]
  
  loadssub <- corr %*% weights_nonorm_sub
  
  # varimax rotation
  vartrans <- varimax(loadssub)
  loadsrot <- loadssub %*% vartrans$rotmat 
  weights_nonorm_sub_rot <- weights_nonorm_sub %*% vartrans$rotmat
  weights_normed_sub_rot <- weights_normed_sub %*% vartrans$rotmat
  
  # normalize rotated weights
  wgtrotnormed <- weights_normed_sub_rot
  for(i in 1:ncomps)
  {
    wgtrotnormed[,i] <- wgtrotnormed[,i]/norm(t(as.matrix(weights_normed_sub_rot[,i])))
  }
  
  # component reliabilites from rotated weights
  relrot <- t(wgtrotnormed) %*% rstar %*% wgtrotnormed %*% 
    solve(t(wgtrotnormed) %*% corr %*% wgtrotnormed)
  relrotvec <- matrix( diag(relrot) )
  
  # change some names to follow SAS macro (may be cleaned up later)
  wgtrotnonorm <- weights_nonorm_sub_rot
  wgtrotnormed <- wgtrotnormed
  
  sink(file=file)
  # printed output
  print('RCA component reliability estimates')
  print(relrot)
  print('')
  
  print('Table of non-normed RCA weights for each component')
  print(wgtrotnonorm)
  print('')
  
  print('Table of normalized RCA weights for each component')
  print(wgtrotnormed)
  print('')
  
  print('Table of orthogonally rotated and normed RCA weights for retained components')
  print(wgtrotnormed)
  print('')
  
  print('Table of orthogonally rotated loadings for retained components')
  print(loadsrot)
  print('')
  
  sink()
  
  # returned output
  list(reliabilities=reliabilities, wgtrotnonorm=wgtrotnonorm,
       wgtrotnormed=wgtrotnormed, loadsrot=loadsrot)
}


