#' Generate random variates from a beta-Poisson distribution
#'
#' @param n The number of random values to return
#' @param alpha A positive value, the alpha shape parameter of the beta distribution
#' estimated by using \code{get_burstiness_bpsc},
#' \code{get_burstiness_scale}, or other method of the user's choosing.
#' @param beta A positive value, the beta shape parameter of the beta distribution
#' estimated by using \code{get_burstiness_bpsc},
#' \code{get_burstiness_scale}, or other method of the user's choosing.
#' @param scale A positive value, the scaling parameter for the beta distribution
#' estimated by using \code{get_burstiness_bpsc},
#' \code{get_burstiness_scale}, or other method of the user's choosing.
#' @return A vector of the random variates
#' @export
#' @examples
#' set.seed(8675309)
#' bp.vec=rBetaPois(100,0.1,1.0,3000)
#' hist(bp.vec, prob=T)
#'
rBetaPois = function(n,alpha,beta,scale=1) {

  if (n <= 0 | !is.numeric(n) | n!=round(n)){
    stop("n must be a positive integer")
  }

  if(alpha<=0 | !is.numeric(alpha)){
    stop("alpha must be a number greater than 0.")
  }

  if(beta<=0 | !is.numeric(beta)){
    stop("beta must be a number greater than 0.")
  }

  if(scale<=0 | !is.numeric(scale)){
    stop("scale must be a number greater than 0.")
  }
  lambda=scale*stats::rbeta(n,shape1=alpha,shape2=beta)
  bp.out=stats::rpois(n,lambda)
  return(bp.out)
}
