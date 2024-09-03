#' Simulate example data for SQL
#'
#' Simulate an additve factor model to illustrate the use of SQL.
#' The generating function are simulated randomly using trigonometric basis functions.
#' The factor are simulated as independent normally distributed variables (unless specified otherwise).
#' 
#' @param n sample size.
#' @param p number of variables.
#' @param q number of factors.
#' @param sde standard deviation of the errors.
#' @param factor a n by q matrix of latent factors. 
#' If not specified, the factors are generated as independent normally distributed variables.
#' @param generator a list of q functions. 
#' Each function should take as input a vector of length m and produce a m by p matrix.
#' If not specified, the factors are generated as using random trigonometric functions.
#' @returns 
#' A list with elements:
#' \itemize{
#' \item \code{data} the data matrix.
#' \item \code{factor} a n times q matrix  containing the latent factors.
#' \item \code{generator} a list of q functions. 
#' Each function takes as input a vector of length m and produce a m by p matrix.
#' \item \code{SNR} Signal to noise ratio.
#' }
#' @export 
#' @examples
#' set.seed(123456)
#' sim <- simulate_afm(n = 150, p = 200, q = 1)
#' sql <- SQL(sim$data )
#' sql
#' plot(sql)
#' abs( cor(sim$factor, sql$factor) )
simulate_afm <- function(n, p, q = 1, sde = 1, factor = NULL, generator = NULL ){
  is_positive_integer <- function(x) (x == as.integer(x)) & (x > 0)
  stopifnot("n should be a positive integer" = is_positive_integer(n) )
  stopifnot("p should be a positive integer" = is_positive_integer(p) )
  stopifnot("q should be a positive integer" = is_positive_integer(q) )
  stopifnot("sde should be positive number" = sde > 0 )
  if(!is.null(factor)){
    stopifnot("factor should be a matrix with q columns" = is.matrix(factor) & ncol(factor) == q )
    stopifnot("factor should have n rows" = nrow(data) == n )
  }
  if(!is.null(generator)){
    stopifnot("generator should be a list of length q" = is.list(generator) & length(generator) == q )
  }
  if(is.null(factor)){
    factor <- matrix( rnorm(n*q ), ncol = q )
    factor <- apply(factor, 2, scale)
    factor <- factor %*% solve( expm::sqrtm( var(factor) ) )    
  }
  if(is.null(generator)){
    generator <- replicate( q,{
      f <- replicate( p , gen_func(rnorm(8)), simplify = FALSE ) 
      function(z){
        matrix( sapply( f, function(g) g(z) ), ncol = p )
      }
    }, simplify = FALSE )
  }
  eps <- matrix( rnorm( n * p, sd = sde ), ncol = p )
  x_expected <- Reduce( '+', lapply( 1:q, function(l) generator[[l]](factor[,l]) ) )
  x <- x_expected + eps
  SNR <- mean( apply( x_expected, 2, var ) / sde^2 )
  return( list( data = x, factor = factor, generator = generator, SNR = SNR ) )
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
  is_positive_integer <- function(x) (x == as.integer(x)) & (x > 0)
  stopifnot("theta should be a vector of even length" = is_positive_integer(length(theta) / 2) )
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
