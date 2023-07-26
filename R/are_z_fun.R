#' Absolute Reconstruction Error for Z
#'
#' Compute the absolute reconstruction error associated with inference of
#' the clonal configuration matrix, Z.
#'
#' @param true.Z an object of class matrix representing the true clonal tree
#' configuration matrix, Z
#' @param inferred.Z an object of class matrix representing the inferred
#' clonal tree configuration matrix, Z
#'
#' @return
#' A numeric representing the absolute reconstruction error associated with
#' inference of the clonal configuration matrix, Z.
#'
#' @details
#' Please refer to Algorithm S3 in the manuscript text for details regarding
#' how the absolute reconstruction error is computed. In brief, the error for
#' Z is quantified by finding the minimum absolute distance between
#' the inference and all permutations of the truth by 1) computing the distance
#' metric, 2) computing extra error from additional columns, and 3) estimating
#' the ARE for each component as a combination of the error from steps 1) and 2).
#'
#' @examples
#' # Simulate alternative and total read counts for bulk and single-cell data
#' sims.out<-simulate_data(N=15, S=2, M=5, alpha=0.1, beta=1, kappa=1, tau=999,
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
#'                         niter=5000, nchains=3, thin=10, pburn=0.1,
#'                         seed=8675309)
#'
#' # Get best tree across all chains and subclones via DIC
#' best.tree.out<-get_best_tree(get.trees.out)
#'
#' # Compute absolute reconstruction error for Z
#' are_z_fun(true.Z=sims.out$true.tree$Z, inferred.Z=best.tree.out$tree$Z)
#'
#' @export

are_z_fun <- function(true.Z, inferred.Z){

  # check arguments
  if (!inherits(true.Z, "matrix")){
    stop("true.Z must be of class \"matrix\"")
  }

  if(!all(true.Z %in% 0:1)){
    stop("true.Z must be a binary matrix")
  }

  if (!inherits(inferred.Z, "matrix")){
    stop("inferred.Z must be of class \"matrix\"")
  }

  if(!all(inferred.Z %in% 0:1)){
    stop("inferred.Z must be a binary matrix")
  }

  if (!(all(rownames(true.Z) == rownames(inferred.Z)))) {
    stop("Rownames for true.Z and inferred.Z are not identical")
  }

  if (!(all(colnames(true.Z) == colnames(inferred.Z)))) {
    stop("Colnames for true.Z and inferred.Z are not identical")
  }

  if (!(identical(dim(true.Z), dim(inferred.Z)))) {
    stop("Dimensions for true.Z and inferred.Z are not identical")
  }

  # Check for duplicated columns in inferred Z and remove if present ONLY
  # if inference is larger than truth (i.e., number of columns in inference
  # exceeds number of columns in truth)
  if(any(duplicated(t(inferred.Z))) & ncol(inferred.Z) > ncol(true.Z)){
    # Which columns are duplicated
    id.Z.dup <- which(duplicated(t(inferred.Z)))
    # Remove duplicates
    inferred.Z <- inferred.Z[,-id.Z.dup]
  }

  M <- nrow(inferred.Z)
  K <- ncol(inferred.Z)

  #--------------------------------
  # Compute distance metric for Z
  #--------------------------------
  dist.Z <- dist_metric_fun(truth=true.Z, inferred=inferred.Z)$dist.metric

  #--------------------------------------------------------------
  # Compute extra error from additional column(s), if applicable
  #--------------------------------------------------------------
  # Extra error for Z
  ncol.diff.Z <- abs(ncol(inferred.Z)-ncol(true.Z))

  # If number of columns in the inference differs from number of columns in
  # the truth
  if(ncol.diff.Z!=0){
    # Add additional error that is equal to the length of the extra
    # column(s), M, times the number of extra column(s)
    e.Z <- M*ncol.diff.Z
    # Otherwise, set error due to extra columns to 0
  }else{e.Z <- 0}

  #-------------------------
  # Estimate the ARE for Z
  #-------------------------

  are.Z <- dist.Z+e.Z

  #--------------------------------------------------
  # Convert the estimated ARE for Z to a percentage
  #--------------------------------------------------

  # Divide by the total number of positions in the largest matrix
  are.Z <- are.Z/max(prod(dim(true.Z)),prod(dim(inferred.Z)))*100

  are.Z
}