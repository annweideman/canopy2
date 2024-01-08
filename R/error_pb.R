#' Absolute Reconstruction Error for Pb
#'
#' Compute the absolute reconstruction error associated with inference of
#' the bulk sample to clone assignment matrix, Pb.
#'
#' @param true.Pb an object of class matrix representing the true bulk sample to
#' clone assignment matrix, Pb
#' @param inferred.Pb an object of class matrix representing the inferred
#' bulk sample to clone assignment matrix, Pb
#'
#' @return
#' A numeric representing the absolute reconstruction error associated with
#' inference of the bulk sample to clone assignment matrix, Pb.
#'
#' @details
#' Please refer to Algorithm S3 in the manuscript text for details regarding
#' how the absolute reconstruction error is computed. In brief, the error for
#' Pb is quantified by finding the minimum absolute distance between
#' the inference and all permutations of the truth by 1) computing the distance
#' metric, 2) computing extra error from additional columns, and 3) estimating
#' the ARE for Pb as a combination of the error from steps 1) and 2).
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
#' # Note: this is written to quickly compile for demonstration purposes, as we
#' # know the optimal number of subclones is 6. In practice, we would attempt a
#' # larger range of subclones (e.g., Klist=3:10) and a larger number of
#' # iterations and chains (e.g., niter=50000 and nchains=10).
#' get.trees.out<-get_trees(Rs=Rs, Rb=Rb, Xs=Xs, Xb=Xb,
#'                         alpha=param.out$alpha, beta=param.out$beta,
#'                         kappa=1, tau=999, Klist=3:5,
#'                         niter=5000, nchains=5, thin=20, pburn=0.2,
#'                         ncores=4, seed=8675309)
#'
#' # Get best tree across all chains and subclones via DIC
#' best.tree.out<-get_best_tree(get.trees.out)
#'
#' # Compute absolute reconstruction error for Ps
#' error_pb(true.Pb=sims.out$true.tree$Pb, inferred.Pb=best.tree.out$tree$Pb)
#'
#' @export

error_pb <- function(true.Pb, inferred.Pb){

  # check arguments
  if (!inherits(true.Pb, "matrix")){
    stop("true.Pb must be of class \"matrix\"")
  }

  if(!all(true.Pb >= 0 & true.Pb <= 1)){
    stop("true.Pb must be a decimal between 0 and 1 inclusive")
  }

  if (!inherits(inferred.Pb, "matrix")){
    stop("inferred.Pb must be of class \"matrix\"")
  }

  if(!all(inferred.Pb >= 0 & inferred.Pb <= 1)){
    stop("inferred.Pb must be a decimal between 0 and 1 inclusive")
  }

  # Transpose to have subclones as columns
  true.Pb <- t(true.Pb); inferred.Pb <- t(inferred.Pb)

  # Cannot evaluate using == 1 due to floating point precision (sums are very
  # slightly different from 1)
  if(!all.equal(as.numeric(rowSums(true.Pb)),rep(1,nrow(true.Pb)))){
    stop("true.Pb must have row sums equal to 1")
  }

  if(!all.equal(as.numeric(rowSums(inferred.Pb)),rep(1,nrow(inferred.Pb)))){
    stop("inferred.Pb must have row sums equal to 1")
  }

  # Compute distance metric
  out.dist.Pb <- dist_metric_fun(truth=true.Pb, inferred=inferred.Pb)
  dist.Pb <- out.dist.Pb$dist.metric
  perm.min.Pb <- out.dist.Pb$mat.perm.min

  # Compute extra error from additional column(s), if applicable
  ncol.diff.Pb <- abs(ncol(inferred.Pb)-ncol(true.Pb))

  # If number of columns in the inference differs from number of columns in
  # the truth
  if(ncol.diff.Pb!=0){

    # Add additional error that is equal to the sum of the elements in the extra
    # column(s)
    if(ncol(inferred.Pb)>ncol(true.Pb)){
      e.Pb <- sum(inferred.Pb[,!(data.frame(inferred.Pb) %in% data.frame(perm.min.Pb))])
    }else{
      e.Pb <- sum(perm.min.Pb[,!(data.frame(perm.min.Pb) %in% data.frame(inferred.Pb))])
    }

    # Otherwise, set error due to extra columns to 0
  }else{e.Pb <- 0}

  # Estimate the ARE for Pb
  are.Pb <- dist.Pb+e.Pb

  # Convert the estimated ARE to a percentage
  are.Pb <- are.Pb/(2*max(nrow(inferred.Pb), nrow(true.Pb)))*100
  return(are.Pb)

}
