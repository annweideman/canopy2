#' Distance Metric
#'
#' Compute the absolute distance between a true matrix and a permuted matrix.
#'
#' @param truth a matrix of binary or continuous (positive) values that
#' represents the ground truth
#' @param inferred a matrix of binary or continuous (positive) values that
#' represents the inference
#'
#' @return
#' A list containing: \code{dist.metric}, a numeric corresponding to the minimum
#' distance metric and \code{mat.perm.min}, a matrix representing the
#' permutation corresponding to the minimum distance metric.
#'
#' @details
#' Please refer to Algorithm S3 in the manuscript text for details regarding
#' how the distance metric is computed. In brief, the distance metric is
#' computed by i) generating all columnwise permutations of whichever matrix
#' (inferred or true) has more columns and ii) for each permutation, computing
#' the distance metric as defined by Algorithm S3.
#'
#' @examples
#' # Simulate read counts
#' sims.out<-simulate_data(N=15, S=2, M=5, alpha=0.1, beta=1, kappa=1, tau=999,
#'                        Ktrue=4, b.mindepth=30, b.maxdepth=50, sc.mindepth=80,
#'                        sc.maxdepth=120, scale=300, seed=8675309)
#'
#' # Run Canopy2 to get list of phylogenetic trees corresponding to all chains
#' # and all subclones
#' get.trees.out<-get_trees(Rs=sims.out$Rs, Rb=sims.out$Rb,
#'                         Xs=sims.out$Xs, Xb=sims.out$Xb,
#'                         alpha=sims.out$alpha, beta=sims.out$beta,
#'                         kappa=1, tau=999, Klist=3:5,
#'                         niter=5000, nchains=3, thin=10, pburn=0.1, seed=8675309)
#'
#' # Get best tree across all chains and subclones via DIC
#' best.tree.out<-get_best_tree(get.trees.out)
#'
#' # Distance between true and inferred clonal configuration matrices, Z
#' dist.Z<-dist_metric_fun(truth=sims.out$true.tree$Z, inferred=best.tree.out$tree$Z)
#'
#' # Print the minimum distance metric
#' dist.Z$dist.metric
#'
#' # Print the permutation corresponding to the minimum distance metric
#' dist.Z$mat.perm.min
#'
#' @export

dist_metric_fun <- function(truth, inferred){

  if (!inherits(truth, "matrix")){
    stop("truth must be of class \"matrix\"")
  }
  if (!inherits(inferred, "matrix")){
    stop("inferred must be of class \"matrix\"")
  }
  if(!is.numeric(truth) | !all(truth>=0)){
    stop("truth must contain non-negative values")
  }
  if(!is.numeric(inferred) | !all(inferred>=0)){
    stop("inferred must contain non-negative values")
  }

  # Find the maximum number of columns
  max.ncol <- max(ncol(truth),ncol(inferred))

  # Store all possible permutations of column numbers
  perms <- gtools::permutations(n=max.ncol,r=max.ncol,v=1:max.ncol)
  # If the matrices differ in size, then remove the overage from the
  # permutations
  perms <- perms[,1:min(ncol(inferred),ncol(truth))]

  # Generate all columnwise permutations of the matrix with the most number
  # of columns
  # If number of columns in the inference exceeds or equals the number of
  # columns in the truth
  if(ncol(truth) <= ncol(inferred)){
    mat.perms <- lapply(1:nrow(perms), function(x) inferred[,perms[x,]])
    # Compute the distance between the truth and each permuted inference
    mat.dist <- lapply(1:length(mat.perms), function(x)
      sum(colSums(abs(truth-mat.perms[[x]]))))
    # If number of columns in the truth exceeds number of columns in the inference
  }else{
    mat.perms <- lapply(1:nrow(perms), function(x) truth[,perms[x,]])
    # Compute the distance between each permuted truth and the inference
    mat.dist <- lapply(1:length(mat.perms), function(x)
      sum(colSums(abs(mat.perms[[x]]-inferred))))
  }

  # Compute the minimum distance metric
  dist.metric <- min(unlist(mat.dist))

  # Store the permutation corresponding to the minimum distance metric
  mat.perm.min <- mat.perms[[which.min(unlist(mat.dist))]]

  return(list(dist.metric=dist.metric,mat.perm.min=mat.perm.min))
}
