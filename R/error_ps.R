#' Absolute Reconstruction Error for Ps
#'
#' Compute the absolute reconstruction error associated with inference of
#' the single-cell to clone assignment matrix, Ps.
#'
#' @param true.Ps an object of class matrix representing the true single-cell to
#' clone assignment matrix, Ps
#' @param inferred.Ps an object of class matrix representing the inferred
#' single-cell to clone assignment matrix, Ps
#'
#' @return
#' A numeric representing the absolute reconstruction error associated with
#' inference of the single-cell to clone assignment matrix, Ps.
#'
#' @details
#' Please refer to Algorithm S3 in the manuscript text for details regarding
#' how the absolute reconstruction error is computed. In brief, the error for
#' Ps is quantified by finding the minimum absolute distance between
#' the inference and all permutations of the truth by 1) computing the distance
#' metric, 2) computing extra error from additional columns, and 3) estimating
#' the ARE for Ps as a combination of the error from steps 1) and 2).
#'
#' @examples
#' # Simulate alternative and total read counts for bulk and single-cell data
#' sims.out<-simulate_data(N=30, S=2, M=5, alpha=0.1, beta=1, kappa=1, tau=999,
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
#' # Run Canopy2 to get list of phylogenetic trees corresponding to all chains
#' # and all subclones
#' get.trees.out<-get_trees(Rs=Rs, Rb=Rb, Xs=Xs, Xb=Xb,
#'                         alpha=param.out$alpha, beta=param.out$beta,
#'                         kappa=1, tau=999, Klist=3:5,
#'                         niter=10000, nchains=20, thin=20, pburn=0.5,
#'                         ncores=4, seed=8675309)
#'
#' # Get best tree across all chains and subclones via DIC
#' best.tree.out<-get_best_tree(get.trees.out)
#'
#' # Compute absolute reconstruction error for Ps
#' error_ps(true.Ps=sims.out$true.tree$Ps, inferred.Ps=best.tree.out$tree$Ps)
#'
#' @export

error_ps <- function(true.Ps, inferred.Ps){

  # check arguments
  if (!inherits(true.Ps, "matrix")){
    stop("true.Ps must be of class \"matrix\"")
  }

  if(!all(true.Ps >= 0 & true.Ps <= 1)){
    stop("true.Ps must be a decimal between 0 and 1 inclusive")
  }

  if (!inherits(inferred.Ps, "matrix")){
    stop("inferred.Ps must be of class \"matrix\"")
  }

  if(!all(inferred.Ps >= 0 & inferred.Ps <= 1)){
    stop("inferred.Ps must be a decimal between 0 and 1 inclusive")
  }

  # Number of single-cells
  N<-ncol(true.Ps)

  # Transpose to have subclones as columns
  true.Ps <- t(true.Ps); inferred.Ps <- t(inferred.Ps)

  # Cannot evaluate using == 1 due to floating point precision (sums are very
  # slightly different from 1)
  if(!all.equal(as.numeric(rowSums(true.Ps)),rep(1,nrow(true.Ps)))){
    stop("true.Ps must have row sums equal to 1")
  }

  if(!all.equal(as.numeric(rowSums(inferred.Ps)),rep(1,nrow(inferred.Ps)))){
    stop("inferred.Ps must have row sums equal to 1")
  }

  # Compute distance metric
  dist.Ps <- dist_metric_fun(truth=true.Ps, inferred=inferred.Ps)$dist.metric

  # Compute extra error from additional column(s), if applicable
  ncol.diff.Ps <- abs(ncol(inferred.Ps)-ncol(true.Ps))

  # If number of columns in the inference differs from number of columns in
  # the truth
  if(ncol.diff.Ps!=0){
    # Add additional error that is equal to the length of the extra
    # column(s), N, times the number of extra column(s)
    e.Ps <- N*ncol.diff.Ps
    # Otherwise, set error due to extra columns to 0
  }else{e.Ps <- 0}

  # Estimate the ARE for Ps
  are.Ps <- dist.Ps+e.Ps

  # Convert the estimated ARE to a percentage
  are.Ps <- are.Ps/max(prod(dim(inferred.Ps)),prod(dim(true.Ps)))*100
  return(are.Ps)

}
