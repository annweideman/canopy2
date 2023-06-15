#' Get clonal composition
#'
#' Returns the single nucleotide variants (SNVs) belonging to each subclone.
#' This function adapts the \code{getclonalcomposition} function from
#' \insertCite{Jiang2017}{Canopy2} by removing any portions of the code related
#' to copy number alterations (CNAs).
#'
#' @param tree phylogenetic tree of class \code{phylo}.
#'
#' @return
#' An M (mutations) x 1 matrix of the mutations assigned to each subclone.
#' The first subclone, "Clone 1" should not have any assigned mutations since it
#' corresponds to the normal subclone.
#'
#' @examples
#' # Load post-processed data for patient GBM10
#' data("GBM10_postproc")
#'
#' # Run Canopy2 to get list of phylogenetic trees corresponding to all chains
#' # and all subclones
#' get.trees.out<-get_trees(Rs=GBM10_postproc@Rs, Rb=GBM10_postproc@Rb,
#'                         Xs=GBM10_postproc@Xs, Xb=GBM10_postproc@Xb,
#'                         alpha=GBM10_postproc@param.est$alpha,
#'                         beta=GBM10_postproc@param.est$beta, kappa=1,
#'                         tau=999, Klist=4:6, niter=1000, nchains=5, thin=10,
#'                         pburn=0.1, seed=8675309)
#'
#' # Examine diagnostic plots
#' get_diagnostics(get.trees.out, project=NULL, outpath=NULL)
#'
#' # Get best tree across all chains and subclones via DIC
#' best.tree.out<-get_best_tree(get.trees.out)
#'
#' # Return the single nucleotide variants (SNVs) belonging to each subclone
#' clonal.muts<-get_clonal_composition(tree=best.tree.out$tree)
#'
#' clonal.muts
#'
#' @export

get_clonal_composition = function(tree) {
  snaname = rownames(tree$snv)
  n = (nrow(tree$edge) + 2)/2
  clonal.muts = vector("list", n)
  for (tip in 2:n) {
    child.node = tip
    parent.node = tree$edge[which(tree$edge[, 2] == child.node), 1]
    while (parent.node >= (n + 1)) {
      muttemp = snaname[intersect(which(tree$snv[, 2] == parent.node),
                                  which(tree$snv[, 3] == child.node))]
      if (length(muttemp) > 0) {
        clonal.muts[[tip]] = c(clonal.muts[[tip]],
                                    muttemp)
      }
      child.node = parent.node
      if (child.node == (n + 1))
        break
      parent.node = tree$edge[which(tree$edge[, 2] == child.node),
                              1]
    }
  }

  clonal.muts[[1]] = "None"
  for (k in 1:length(clonal.muts)) {
    if (!is.null(clonal.muts[[k]]) & (length(clonal.muts[[k]])) >
        0) {
      clonal.muts[[k]] = sort(clonal.muts[[k]])
    }
  }

  clonal.muts<-lapply(clonal.muts,paste,collapse=", ")
  names(clonal.muts)<-paste0("Clone:",1:length(clonal.muts))

  clonal.muts<-do.call('rbind',clonal.muts)
  return(clonal.muts)
}
