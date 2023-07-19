#' The Beta-Poisson Distribution
#'
#' This function is used to generate random variants for the beta-Poisson
#' distribution. The beta distribution is modified to have a support of
#' \[0, \code{scale}\] instead of the typical \[0, 1\], and \code{alpha} and
#' \code{beta} represent the shape parameters.
#'
#' @param n the number of random values to return
#' @param alpha a positive value, the alpha shape parameter of the beta distribution
#' estimated by using \code{get_burstiness_bpsc},
#' \code{get_burstiness_scale}, or other method of the user's choosing
#' @param beta a positive value, the beta shape parameter of the beta distribution
#' estimated by using \code{get_burstiness_bpsc},
#' \code{get_burstiness_scale}, or other method of the user's choosing.
#' @param scale a positive value, the scaling parameter for the beta distribution
#' estimated by using \code{get_burstiness_bpsc},
#' \code{get_burstiness_scale}, or other method of the user's choosing
#'
#' @return a vector of the random variates from the beta-Poisson distribution
#'
#' @details
#' Each of the beta-Poisson random variates, \eqn{r}, are generated as
#'
#' \deqn{r \sim Pois(\lambda),}
#'
#' where \eqn{\lambda \sim c*Beta(\alpha,\beta)} for \eqn{c=}\code{scale},
#' \eqn{\alpha}=\code{alpha}, and \eqn{\beta}=\code{beta}.
#'
#' @examples
#' set.seed(8675309)
#' bp.vec=rBetaPois(100,0.1,1.0,3000)
#' hist(bp.vec, prob=TRUE)
#'
#' @export

rBetaPois = function(n,alpha,beta,scale=1) {

  if (n <= 0 | !is.numeric(n) | n!=round(n)){
    stop("n must be a positive integer")
  }

  if(alpha<=0 | !is.numeric(alpha)){c
    stop("alpha must be a number greater than 0.")
  }

  if(beta<=0 | !is.numeric(beta)){
    stop("beta must be a number greater than 0.")
  }

  if(scale<=0 | !is.numeric(scale)){
    stop("scale must be a number greater than 0.")
  }
  lambda<-scale*stats::rbeta(n,shape1=alpha,shape2=beta)
  stats::rpois(n,lambda)
}
