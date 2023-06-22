#' Absolute Reconstruction Error for Z, Ps, and Pb
#'
#' Compute the absolute reconstruction error (ARE) associated with inference of
#' the clonal configuration matrix, Z, cell-to-clone assignment matrix, Ps,
#' and sample-to-clone assignment matrix, Pb.
#'
#' @param true.tree an object of class \code{phylo} representing the ground truth
#' phylogenetic tree
#' @param inferred.tree an object of class \code{phylo} representing the inferred
#' phylogenetic tree
#'
#' @return
#' A list containing: \code{are.Z}, the absolute reconstruction error associated
#' with inference of the clonal configuration matrix, Z, \code{are.Ps}, the
#' absolute reconstruction error associated with inference of the cell-to-clone
#' assignment matrix, Ps, and \code{are.Pb}, the absolute reconstruction error
#' associated with inference of the sample-to-clone assignment matrix, Pb.
#'
#' @details
#' Please refer to Algorithm S3 in the manuscript text for details regarding
#' how the absolute reconstruction error is computed. In brief, the error for
#' each component is quantified by finding the minimum absolute distance between
#' the inference and all permutations of the truth by 1) computing the distance
#' metric for each component, 2) computing extra error from additional columns,
#' and 3) estimating the ARE for each component as a combination of the error
#' from steps 1) and 2).
#'
#' @examples
#' # Simulate alternative and total read counts for bulk and single-cell data
#' sims.out<-simulate_data(N=15, S=2, M=5, alpha=0.1, beta=1, kappa=1, tau=999,
#'                        Ktrue=4, b.mindepth=30, b.maxdepth=50, sc.mindepth=80,
#'                        sc.maxdepth=120, scale=300, seed=8675309)
#'
#' # Run Canopy2 to get list of phylogenetic trees corresponding to all chains
#' # and all subclones
#' get.trees.out<-get_trees(Rs=sims.out$Rs, Rb=sims.out$Rb,
#'                          Xs=sims.out$Xs, Xb=sims.out$Xb,
#'                          alpha=sims.out$alpha, beta=sims.out$beta,
#'                          kappa=1, tau=999, Klist=3:5,
#'                          niter=5000, nchains=3, thin=10, pburn=0.1,
#'                          seed=8675309)
#'
#' # Get best tree across all chains and subclones via DIC
#' best.tree.out<-get_best_tree(get.trees.out)
#'
#' # Compute absolute reconstruction error for each component
#' are.out<-are_fun(true.tree=sims.out$true.tree, inferred.tree=best.tree.out$tree)
#'
#' # Print absolute reconstruction error for each component
#' are.out$are.Z; are.out$are.Ps; are.out$are.Pb
#'
#' @export

are_fun <- function(true.tree, inferred.tree){

  # check arguments
  if (!inherits(true.tree, "phylo")){
    stop("true.tree must be of class \"phylo\"")
  }
  if (!inherits(inferred.tree, "phylo")){
    stop("inferred.tree must be of class \"phylo\"")
  }
  if (is.null(true.tree$Z)){
    stop("true.tree must contain clonal configuration matrix, Z, under
         slot true.tree$Z")
  }
  if (is.null(true.tree$Ps)){
    stop("true.tree must contain cell-to-clone assignment matrix, Ps, under
         slot true.tree$Ps")
  }
  if (is.null(true.tree$Pb)){
    stop("true.tree must contain sample-to-clone assignment matrix, Pb, under
         slot true.tree$Pb")
  }
  if (is.null(inferred.tree$Z)){
    stop("inferred.tree must contain clonal configuration matrix, Z, under
         slot inferred.tree$Z")
  }
  if (is.null(inferred.tree$Ps)){
    stop("inferred.tree must contain cell-to-clone assignment matrix, Ps, under
         slot inferred.tree$Ps")
  }
  if (is.null(inferred.tree$Pb)){
    stop("inferred.tree must contain sample-to-clone assignment matrix, Pb, under
         slot inferred.tree$Pb")
  }

  N<-ncol(inferred.tree$Ps)
  M<-nrow(inferred.tree$Z)

  #--------------------------------------------------------------
  # Transpose when needed to have subclones as columns
  #--------------------------------------------------------------
  true.Z <- true.tree$Z; inferred.Z <- inferred.tree$Z
  true.Ps <- t(true.tree$Ps); inferred.Ps <- t(inferred.tree$Ps)
  true.Pb <- t(true.tree$Pb); inferred.Pb <- t(inferred.tree$Pb)

  #-------------------------------------------
  # Compute distance metric for each component
  #-------------------------------------------
  dist.Z <- dist_metric_fun(truth=true.Z, inferred=inferred.Z)$dist.metric
  dist.Ps <- dist_metric_fun(truth=true.Ps, inferred=inferred.Ps)$dist.metric
  out.dist.Pb <- dist_metric_fun(truth=true.Pb, inferred=inferred.Pb)
  dist.Pb <- out.dist.Pb$dist.metric
  perm.min.Pb <- out.dist.Pb$mat.perm.min

  #--------------------------------------------------------------
  # Compute extra error from additional column(s), if applicable
  #--------------------------------------------------------------
  ### Extra error for Z
  ncol.diff.Z <- abs(ncol(inferred.Z)-ncol(true.Z))

  # If number of columns in the inference differs from number of columns in
  # the truth
  if(ncol.diff.Z!=0){
    # Add additional error that is equal to the length of the extra
    # column(s), M, times the number of extra column(s)
    e.Z <- M*ncol.diff.Z
    # Otherwise, set error due to extra columns to 0
  }else{e.Z <- 0}

  ### Extra error for Ps
  ncol.diff.Ps <- abs(ncol(inferred.Ps)-ncol(true.Ps))

  # If number of columns in the inference differs from number of columns in
  # the truth
  if(ncol.diff.Ps!=0){
    # Add additional error that is equal to the length of the extra
    # column(s), N, times the number of extra column(s)
    e.Ps <- N*ncol.diff.Ps
    # Otherwise, set error due to extra columns to 0
  }else{e.Ps <- 0}

  ### Extra error for Pb
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

  #-------------------------------------------
  # Estimate the ARE for each component
  #-------------------------------------------

  are.Z <- dist.Z+e.Z
  are.Ps <- dist.Ps+e.Ps
  are.Pb <- dist.Pb+e.Pb

  #-------------------------------------------
  # Convert the estimated AREs to percentages
  #-------------------------------------------

  are.Z <- are.Z/max(prod(dim(inferred.Z)),prod(dim(true.Z)))*100
  are.Ps <- are.Ps/max(prod(dim(inferred.Ps)),prod(dim(true.Ps)))*100
  are.Pb <- are.Pb/(2*max(nrow(inferred.Pb), nrow(true.Pb)))*100

  return(list(are.Z=are.Z, are.Ps=are.Ps, are.Pb=are.Pb))

}
