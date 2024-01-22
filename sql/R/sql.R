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
#' @param use_Rfast If TRUE, use package Rfast to compute the cross-product of the data matrix. This option may lead to speed gains if the number of variables is large.
#' @param pretrain_tol Tolerance for pretraining.
#' @param Sigma Cross-product of the centered data.
#' @returns 
#' An object with S3 class sql.
#' * factor a n times q matrix  containing the estimated latent factors.
#' * data the data matrix.
#' * qap a list containing the solution of the Quadratic Assignment Problem, the mean squared error for each iteration, and the convergence status.
#' * splines an object to produce the basis functions.
#' * lambda the tuning parameter.
#' * d dimension of the basis functions.
#' * q number of latent factors.
#' * intercepts parameter vector of intercepts of length p.
#' * EV Explained variance.
#' @export 
#' @examples
#' set.seed(123456)
#' sim <- simulate_afm(n = 150, p = 200, q = 1, sde = 1)
#' sql <- SQL(sim$data, d = 4 )
#' sql
#' plot(sql)
#' abs( cor(sim$factor_z, sql$factor) )
SQL <- function( data, q = 1, d = 4, lambda = 0.1, tol = 1e-4, max_iter = 30, max_cycles = 10, use_Rfast = FALSE, 
                 pretrain_tol = 1e-4, Sigma = NULL ){
  # check args:
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
  if(!is.null(Sigma)){
    stopifnot("Sigma should be a diagonal positive definite matrix" = is.positive.definite(Sigma) )
    stopifnot(nrow(Sigma) == nrow(data) )
    if(q>1) warning("Argument Sigma is ignored")
  }
  if( ncol(data) < 4 ) warning("number of columns might be too small")
  if( ncol(data) < q * d ) warning("for identifiability q times d shoul be less than number of columns.")
  # parameter list:
  par <- list( n= nrow(data), p = ncol(data), lambda = lambda, tol = tol, max_iter = max_iter, 
               max_cycles = max_cycles, use_Rfast = use_Rfast, pretrain_tol = pretrain_tol,
               print = TRUE)
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


#' Simulate example data for Additive Factor Model
#'
#' Function used to simulated data sets to illustrate the use of SQL.
#' The generating function are simulated randomly using trigonometric basis functions.
#' The factor are simulated as independent normally distributed variables (unless specified otherwise).
#' 
#' @param n sample size.
#' @param p number of variables.
#' @param q number of factors.
#' @param sde standard deviation of the errors.
#' @param factor a n by q matrix of latent factors. 
#' If not specified, the factors are generated as independent normally distributed variables.
#' @param generator a list of length q functions. 
#' Each function takes as input a vector of length m and produce a m by p matrix.
#' If not specified, the factors are generated as using random trigonometric functions.
#' @export
simulate_afm <- function(n, p, q = 1, sde = 1, factor = NULL, generator = NULL ){
  # generate factors, perfectly independent with mean 0 and not correlated
  if(is.null(factor)){
    factor <- matrix( rnorm(n*q ), ncol = q )
    factor <- apply(factor, 2, scale)
    factor <- factor %*% solve( expm::sqrtm( var(factor) ) )    
  }
  # functions:
  if(is.null(generator)){
    generator <- replicate( q,{
      f <- replicate( p , gen_func(rnorm(8)), simplify = FALSE ) 
      function(z){
        matrix( sapply( f, function(g) g(z) ), ncol = p )
      }
    }, simplify = FALSE )
  }
  # model:
  eps <- matrix( rnorm( n * p, sd = sde ), ncol = p )
  x_expected <- Reduce( '+', lapply( 1:q, function(l) generator[[l]](factor[,l]) ) )
  x <- x_expected + eps
  SNR <- mean( apply( x_expected, 2, var ) / sde^2 )
  return( list( data = x, factor_z = factor, generator = generator, SNR = SNR ) )
}


#' Simulate example data for Additive Factor Model
#'
#' Generate a function from linear combinations of trigonometric functions.
#' The functions are centered, i.e. E(g(Z)) = 0, if Z is normally distributed.
#' 
#' @param theta vector of coefficients for the basis functions. Length of the vector should be even.
#' @export
#' @examples 
#' g <- gen_func(rnorm(8))
#' z <- seq(-3, 3, l = 100)
#' plot(g(z) ~ z, type = "l")
gen_func <- function(theta){
  theta <- theta / mean( theta^2 )
  m <- length(theta) / 2
  ind <- 1:m
  a <- theta[ind] / ind 
  b <- theta[ m + ind] / ind 
  f <- Vectorize( function(z){ sum( a * cos( 2 * pi * ind * z / 8 ) + b * sin( 2 * pi * ind * z / 8 ) ) / 2  } )
  grid <- qnorm(1:200/201)
  mean_f <- mean( f(grid) )
  return( function(z) f(z) - mean_f )
}



#===================
# High level q = 1:
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
  opt <- qap_solve(Sigma, M, P[[1]], par$p, par$tol, par$max_iter, par$print )
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
    qap <- qap_solve(Sigma, M, P[[l]], par$p, par$tol, par$max_iter, par$print )
    P[[l]] <- qap$P
  }
  return(P)
}


#==============
# Low level:
#==============

qap_solve <- function( Sigma, M, P, p, tol = 1e-4, max_iter = 50, print = TRUE ){
  mse <- list()
  convergence <- FALSE
  for( it in 1:max_iter ){
    P <- hungarian_update(Sigma, M, P)
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
  xpred <- Reduce("+", lapply(1:q, function(l) P[[l]] %*% G[[l]]) )
  mean( ( x - xpred )^2 )
}

get_mse_1 <- function( Sigma, M, Pmat, p){
  # Pmat should be a matrix
  PM <- as.matrix( Pmat %*% t(M) )
  G2 <- t(PM) %*% Sigma %*% PM
  out <- Sigma + G2 - 2 * as.matrix( Pmat ) %*% t(PM) %*% Sigma
  return( mean(diag(out) )/ p )
}

hungarian_update <- function( Sigma, M, Pmat ){
  PM <- as.matrix(Pmat %*% t(M))
  G2 <- t(PM) %*% Sigma %*% PM
  cost <- - 2 * Sigma %*% PM
  cost <- t( t( cost + diag(Sigma) ) + diag(G2) )
  ord <- clue::solve_LSAP( cost )
  return( Matrix::t(get_permutationMatrix( ord ) ) )
}

pred_G <- function (x, P, d, lambda ){
  q <- length(P)
  grid <- 1:nrow(x) / (nrow(x)+1)
  splines <- get_splines( d = d + 1 )
  basis <- get_basis(splines, grid, center = TRUE )
  psi <- basis$psi
  Omega_tot <- as.matrix( do.call(Matrix::bdiag, replicate(q, basis$Omega, simplify = F ) ) )
  psimat <- as.matrix( do.call(cbind, lapply(P, function(A) A %*% psi)) )
  Gram_matrix <- t(psimat) %*% psimat + lambda * Omega_tot
  is_singular <- matrixcalc::is.singular.matrix(as.matrix(Gram_matrix))
  if(is_singular){
    print("Matrix is singular!!")
    Gram_matrix <- Gram_matrix + (lambda/2) * diag(q*d)
  }
  B <- solve(Gram_matrix, t(psimat)) %*% x
  G <- lapply(1:q, function(l) psi %*% B[(l - 1) * d + 1:d, ])
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


