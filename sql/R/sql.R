#' Statistical Quantile Learning
#'
#' Estimate the Additive Factor Model with the Statistical Quantile Learning method.
#' The implementation makes use of penalyzed B-splines basis functions.
#' 
#' @param data  A  real data matrix.
#' @param q  The number of factors to be estimated.
#' @param d The number of B-splines basis functions.
#' @param lambda Tuning constant for penalyzation.
#' @param tol Tolerance for the qap optimization problem.
#' @param max_iter  Maximum number of iterations for the qap algorithm.
#' @param max_cycles  Maximum number for the cycles for the backfitting algorithm.
#' @param use_Rfast If TRUE, use package Rfast to compute the cross-product of the data matrix. This option may lead to speed gains if the number of variables is large. To use it, make sure that Rfast is installed.
#' @param pretrain_tol Tolerance for pretraining.
#' @param Sigma Cross-product of the centered data.
#' @param greedy If TRUE, use a greedy algorithm instead of the Hungarian algorithm. This option leads to large speed gains if the number of observations is large.
#' @returns 
#' A list with class \code{sql}.
#' \itemize{
#' \item \code{factor} a n times q matrix  containing the estimated latent factors.
#' \item \code{data} the data matrix.
#' \item \code{opt} a list containing the solution of the Quadratic Assignment Problem, the mean squared error for each iteration, and the convergence status.
#' \item \code{lambda} the tuning parameter.
#' \item \code{d} dimension of the basis functions.
#' \item \code{q} number of latent factors.
#' \item \code{intercepts} parameter vector of intercepts of length p.
#' \item \code{EV} Explained variance.
#' }
#' @export 
#' @import matrixcalc
#' @import Matrix
#' @examples
#' set.seed(123456)
#' sim <- simulate_afm(n = 150, p = 200)
#' sql <- SQL(sim$data )
#' sql
#' plot(sql)
#' abs( cor(sim$factor, sql$factor) )
#' #==============
#' # q= 3 factor:
#' #==============
#' q <- 3
#' sim <- simulate_afm(n = 150, p = 200, q = q)
#' sql <- SQL(sim$data, q = q)
#' sql
#' abs( cor(sim$factor_z, sql$factor) )
SQL <- function( data, q = 1, d = 4, lambda = 0.1, tol = 1e-4, max_iter = 30, max_cycles = 20, use_Rfast = FALSE, 
                 pretrain_tol = 5e-4, Sigma = NULL, greedy = FALSE){
  # Check args:
  is_positive_integer <- function(x) (x == as.integer(x)) & (x > 0)
  stopifnot("data should be a matrix" = is.matrix(data) & ncol(data) > 1 & nrow(data) > 1)
  stopifnot("q should be a positive integer" = is_positive_integer(q) )
  stopifnot( "d should be an integeger greater or equal to 4" = (d >= 4 & is_positive_integer(d) ) )
  stopifnot(lambda >= 0)
  stopifnot(tol >= 0)
  stopifnot(is_positive_integer(max_iter) )
  stopifnot(is_positive_integer(max_cycles))
  stopifnot(is.logical( use_Rfast ) )
  stopifnot(pretrain_tol>=0)
  stopifnot(is.logical(greedy))
  if(use_Rfast){
    if (!requireNamespace("Rfast", quietly = TRUE)) {
      warning("The Rfast package must be installed. Ignoring use_Rfast argument.")
      use_Rfast <- FALSE
    }
  }
  if(!is.null(Sigma)){
    stopifnot("Sigma should be a diagonal positive definite matrix" = is.positive.definite(Sigma) )
    stopifnot(nrow(Sigma) == nrow(data) )
    if(q>1) warning("Argument Sigma is ignored")
  }
  if( ncol(data) < 4 ) warning("number of columns might be too small")
  if( ncol(data) < q * d ) warning("for identifiability q times d shoul be less than number of columns.")
  # Parameter list:
  par <- list( n= nrow(data), p = ncol(data), lambda = lambda, tol = tol, max_iter = max_iter, 
               max_cycles = max_cycles, use_Rfast = use_Rfast, pretrain_tol = pretrain_tol,
               print = TRUE, greedy = greedy)
  # Initialize:
  intercepts <- colMeans(data)
  x <- scale( data, scale = FALSE)
  P <- get_initial_P(x, q )
  if(q == 1 & is.null(Sigma)) Sigma <- crossprod2(x, par$use_Rfast)
  # Algorithm:
  if(par$pretrain_tol != 0 & d >4 ){
    init <- pretrain(x, P, d, Sigma, par)
    P <- init$P
  }else{init <- NULL}
  opt <- SQL_solve(x, P, d, Sigma, par)
  # Output:
  mse <- tail( opt$mse, 1 )
  EV <- 1 - mse / mean(x^2)
  grid <- 1:nrow(x) / (nrow(x)+1)
  factor <- sapply( opt$P, function(x) qnorm( as.matrix(x %*% grid) ) )
  colnames(factor) <- paste0("Z", 1:q)
  output <- list( factor = factor, data = data, opt = opt, lambda = lambda, d = d, q = q, intercepts = intercepts, EV = EV, 
                  init = init )
  class(output) <- "sql"
  return(output)
}



#===================
# High level:
#===================

SQL_solve <- function(x, P, d, Sigma, par){
  q <- length(P)
  if(q>1){
    opt <- backfitting(x, P, d, par )
  }else{
    opt <- SQL_q1(Sigma, P, d, par )
  }
  return(opt)
}

pretrain <- function(x, P, d, Sigma, par){
  d_index <- 4:d
  par$tol <- par$pretrain_tol
  mse <- list()
  for( i in 1:length(d_index) ){
    opt <- SQL_solve(x, P, d = d_index[i], Sigma, par )
    P <- opt$P
    mse[[i]] <- opt$mse
  }
  return( list(P = P, mse = mse ) )
}


SQL_q1 <- function(Sigma, P, d, par ){
  grid <- 1:par$n / (par$n+1)
  splines <- get_splines( d = d + 1 )
  basis <- get_basis(splines, grid, center = TRUE )
  M <- get_projection_matrix( basis$psi, basis$Omega, par$lambda )
  opt <- qap_solve(Sigma, M, P[[1]], par$p, par$tol, par$max_iter, par$print, par$greedy )
  opt$P <- list( opt$P )
  return(opt)
}


backfitting <- function(x, P, d, par){
  mse <- list()
  convergence <- FALSE
  for(it in 1:par$max_cycles){
    G <- pred_G(x, P, d, par$lambda)
    P <- backfitting_cycles(x, G, P, d, par)
    mse[[it]] <- get_mse(x, G, P)
    if(it>1){
      convergence <- abs( mse[[it]] - mse[[it-1]] ) / mse[[it-1]] < par$tol
      if(convergence) break
    }
  }
  if(!convergence) print('Maximum number of cycles reached without convergence')
  return( list(P = P, mse = unlist( mse ), convergence = convergence ) )
}

backfitting_cycles <- function(x, G, P, d, par ){
  q <- length(P)
  grid <- 1:par$n / (par$n+1)
  splines <- get_splines( d = d + 1 )
  basis <- get_basis(splines, grid, center = TRUE )
  M <- get_projection_matrix( basis$psi, basis$Omega, par$lambda )
  for(l in 1:q ){
    U <- as.matrix( x - Reduce( '+', lapply( (1:q)[-l], function(k) P[[k]] %*% G[[k]] ) ) )
    Sigma <- crossprod2(U, par$use_Rfast)
    qap <- qap_solve(Sigma, M, P[[l]], par$p, par$tol, par$max_iter, par$print, par$greedy )
    P[[l]] <- qap$P
  }
  return(P)
}


#==============
# Low level:
#==============

qap_solve <- function( Sigma, M, P, p, tol = 1e-4, max_iter = 50, print = TRUE, greedy = FALSE){
  mse <- list()
  convergence <- FALSE
  for( it in 1:max_iter ){
    P <- hungarian_update(Sigma, M, P, greedy)
    mse[[it]] <- get_mse_1( Sigma, M, P, p )
    if(it > 1){
      convergence <- ( abs( mse[[it]] - mse[[it-1]] ) / mse[[it-1]] ) < tol
      if( convergence ) break
    }
  }
  if(!convergence) if(print) print('Maximum number of iterations reached without convergence')
  out <- list(P = P, mse = unlist( mse ), convergence = convergence )
  return(out)
}


crossprod2 <- function(x, use_Rfast ){
  if( use_Rfast ){
    return( Rfast::Tcrossprod(x, x) )
  }else{
    return( x %*% t(x) )
  }  
}


get_mse <- function (x, G, P){
  # P is a list of matrices
  q <- length(P)
  xpred <- as.matrix( Reduce("+", lapply(1:q, function(l) P[[l]] %*% G[[l]]) ) )
  mean( ( x - xpred )^2 )
}

get_mse_1 <- function( Sigma, M, Pmat, p){
  # Pmat should be a matrix
  PM <- as.matrix( Pmat %*% t(M) )
  G2 <- t(PM) %*% Sigma %*% PM
  out <- Sigma + G2 - 2 * as.matrix( Pmat ) %*% t(PM) %*% Sigma
  return( mean(diag(out) )/ p )
}

hungarian_update <- function( Sigma, M, Pmat, greedy ){
  PM <- as.matrix(Pmat %*% t(M))
  G2 <- t(PM) %*% Sigma %*% PM
  cost <- - 2 * Sigma %*% PM
  cost <- t( t( cost + diag(Sigma) ) + diag(G2) )
  ord <- solve_LSAP2(cost, greedy)
  return( Matrix::t(get_permutationMatrix( ord ) ) )
}

solve_LSAP2 <- function(cost, greedy = FALSE){
  if(greedy){
    ord <- solve_LSAP_greedy(cost)
  }else{
    ord <- clue::solve_LSAP( cost )
  }
  return(ord)
}

get_Beta <- function(x, P, lambda, d){
  q <- length(P)
  grid <- 1:nrow(x) / (nrow(x)+1)
  basis <- get_basis( get_splines(d+1), grid, center = TRUE )
  Omega_tot <- as.matrix( do.call(Matrix::bdiag, replicate(q, basis$Omega, simplify = F ) ) )
  psimat <- as.matrix( do.call(cbind, lapply(P, function(A) A %*% basis$psi ) ) )
  Gram_matrix <- t(psimat) %*% psimat + lambda * Omega_tot
  is_singular <- matrixcalc::is.singular.matrix(as.matrix(Gram_matrix))
  if(is_singular){
    print("Matrix is singular!!")
    Gram_matrix <- Gram_matrix + (lambda/2) * diag(q*d)
  }
  return( solve(Gram_matrix, t(psimat)) %*% x )
}


pred_G <- function (x, P, d, lambda ){
  # uniform scale
  q <- length(P)
  grid <- 1:nrow(x) / (nrow(x)+1)
  basis <- get_basis(get_splines(d+1), grid, center = TRUE )
  B <- get_Beta(x, P, lambda, d)
  G <- lapply(1:q, function(l) basis$psi %*% B[(l - 1) * d + 1:d, ])
  return(G)
}


get_projection_matrix <- function( psi, Omega, lambda ){
  Gram_matrix <- t( psi ) %*% psi + lambda * Omega
  return( psi %*% solve( Gram_matrix, t( psi ) ) )
}


get_initial_P <- function( x, q ){
  uniformize <- function( x ){ rank( x, ties.method = "random" ) / ( length(x) + 1 ) }
  pca <- irlba::prcomp_irlba( x, n = q, scale = F, center = F )
  factor <- apply( as.matrix( pca$x[, 1:q]), 2, uniformize ) 
  ord <- apply( factor, 2, order )
  P <- apply( ord, 2, get_permutationMatrix )
}


get_permutationMatrix <- function( ord ) Matrix::t( Matrix::sparseMatrix( seq_along( ord ), ord, x = 1 ) )


get_splines <- function( d, range = c(0,1), degree = 3 ){
  nknots <- d - 2
  knots <- orthogonalsplinebasis::expand.knots( seq( range[1], range[2], length.out = nknots ) )
  basis  <-  orthogonalsplinebasis::SplineBasis(knots = knots, order = degree + 1 )
  d_basis <- orthogonalsplinebasis::deriv( basis)
  d2psi <- orthogonalsplinebasis::OuterProdSecondDerivative(basis)
  return( list( basis = basis, d_basis = d_basis, d2psi = d2psi ) )
}

get_basis <- function(splines, grid, center = TRUE ){
  # when centered splines are used, we remove the last basis to ensure linear independence
  psi <- orthogonalsplinebasis::evaluate( splines$basis, grid )
  if(center){
    psi <- scale( psi, scale = FALSE )
    psi <- psi[, 1:(ncol(psi) - 1 ) ]
  }
  d <- ncol(psi)
  Omega <- splines$d2psi[1:d, 1:d]
  return(list(psi = psi, Omega = Omega) )
}



solve_LSAP_greedy <- function(cost){
  rank( apply( cost, 1, which.min ), ties.method = "random" )
}



#==============
# SQL class :
#==============


#' Print a Statistical Quantile Learning
#'
#' The default print method for a sql object
#' 
#' @param object  a fitted sql object produced by SQL().
#' @export
print.sql <- function(object){
  prt <- paste0("SQL object \n q = ", object$q, " factor(s) \n lambda = ", object$lambda, 
                "\n Explained variance = ", round( object$EV, 3), "\n Convergence = ", object$opt$convergence )
  cat( prt )
}




#' Plot for Statistical Quantile Learning
#'
#' Takes a fitted sql object produced by SQL() and plot the factors.
#' 
#' @param object  a fitted sql object produced by SQL().
#' @param ... other graphics parameters to pass on to the plotting commands.
#' If this is not provided, predictions are perform over a grid of m=n normal quantiles.
#' @export
plot.sql <- function(object, ...){
  q <- object$q
  par(mfrow = c(q, 1))
  for(l in 1:q){
    plot(object$factor[, l], type = "l", main = paste( "Factor" , l ), ylab = "", ... ) 
  }
}



#' Predictions for Statistical Quantile Learning
#'
#' Obtain predictions from a fitted sql object produced by SQL().
#' 
#' @param object  a fitted sql object produced by SQL().
#' @param new_data a matrix, with q columns, containing the values (of the latent variables) at which predictions are required.
#' @param type type of predictions, with choices "terms" or "response".
#' @param scaling A probability distribution function corresponding to the distribution of new_data.
#' @export
predict.sql <- function( object, new_data, type = c("terms", "response"), scaling = pnorm ){
  stopifnot("new_data should be a matrix with q columns" = is.matrix(new_data) & ncol(new_data) == object$q )
  type <- match.arg(type)
  # For now it is only available for the Gaussian scale
  d <- object$d
  new_grid <- apply(new_data, 2, scaling )
  splines <- sql:::get_splines(d+1)
  B <- sql:::get_Beta(object$data, object$opt$P, object$lambda, d)
  predicted_terms <- lapply(1:object$q, function(l){
    basis <- sql:::get_basis(splines, new_grid[,l], center = TRUE )
    return( basis$psi %*% B[(l - 1) * d + 1:d, ] )
  })
  if (type == "terms") {
    return(predicted_terms)    
  } else if (type == "response") {
    out <- Reduce('+', predicted_terms)
    return(out)
  }
}
