#' Log beta-binomial density
#'
#' Compute the natural log-transformed beta-binomial density.
#'
#' @param r a matrix of non-negative integers representing the alternative read
#' counts.
#' @param x a matrix of non-negative integers representing the total number of
#' read counts
#' @param shape1 a non-negative numeric vector corresponding to the \eqn{\alpha}
#' shape parameter in the beta-binomial distribution. Must have length equal
#' to the number of rows (number of mutations) in \code{r}.
#' @param shape2 a non-negative numeric vector corresponding to the \eqn{\beta}
#' shape parameter in the beta-binomial distribution. Must have length equal to
#' the number of rows (number of mutations) in \code{r}.
#'
#' @return
#' A numeric representing the log beta-binomial density.
#'
#' @details
#' The beta-binomial distribution has density
#' \deqn{f(r) = \binom{x}{r} \frac{\mathrm{B}(r+\alpha,x-r+\beta)}
#' {\mathrm{B}(\alpha,\beta)},} where \eqn{\Beta(x,y)=\frac{\Gamma(x)\,
#' \Gamma(y)}{\Gamma(x+y)}} is the beta function. The final, log-transformed
#' density is written as
#'
#' \deqn{\ell(r) = log(x!) - log[(x-r)!] - log(r!) + log(|\mathrm{B}(r+\alpha,x-r+\beta)|) -
#' log(|\mathrm{B}(\alpha,\beta)|),}
#'
#' where the log transformation is applied to the absolute value of each beta
#' function since the log of a negative value is undefined. Note that the beta
#' function is only defined in \eqn{\mathbb{R}} for non-negative \eqn{\alpha} and
#' \eqn{\beta} and is infinite if either is zero.
#'
#' @examples
#' # Simulate read counts
#' sims.out<-simulate_data(N=15, S=2, M=5, alpha=0.1, beta=1.0, kappa=1, tau=999,
#'                         Ktrue=4, b.mindepth=30, b.maxdepth=50, sc.mindepth=80,
#'                         sc.maxdepth=120, scale=300, seed=8675309)
#'
#' # Estimate parameters for gene kinetics using the BPSC methodology
#' # Note: can also use get_burstiness_scale() for large datasets, but estimates
#' # are not as reliable
#' param.out<-get_burstiness_bpsc(counts=sims.out$G)
#'
#' # Subset the read counts to include only those mutations with estimable gene
#' # kinetics (all are estimable here, but included for utility when using real
#' # data)
#' Rs<-sims.out$Rs[param.out$id.g,] # single cell alternative read counts
#' Xs<-sims.out$Xs[param.out$id.g,] # single cell total read counts
#' Rb<-sims.out$Rb[param.out$id.g,] # bulk alternative read counts
#' Xb<-sims.out$Xb[param.out$id.g,] # bulk total read counts
#'
#' # Compute log beta-binomial density
#' logdBetaBinom(r=Rs, x=Xs, shape1=param.out$alpha, shape2=param.out$beta)
#'
#' @keywords internal

logdBetaBinom=function(r, x, shape1, shape2){
  # Compute density
  lchoose(x,r)+lbeta(r+shape1, x-r+shape2)-lbeta(shape1, shape2)
}

