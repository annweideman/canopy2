#' Plot Z
#'
#' Plot the clonal tree configuration matrix, Z.
#'
#' @param tree an object of class \code{pyhlo} representing the phylogenetic
#' tree.
#'
#' @return
#' A plot of class \code{ggtree} that is automatically plotted to the graphics
#' device.
#'
#' @details
#' The plotted tree is arranged upright such that time runs vertically from
#' the root to the tips of the tree. The leftmost branch corresponds to the
#' population of normal cells, while the remaining branches correspond to the
#' mutant cells. Each red text box, other than the root, denotes a cluster of
#' point mutations (single-nucleotide variants). For example, "M1: 2 SNVs" would
#' imply that cluster M1 contains two point mutations, perhaps mutations in genes
#' ABC and XYZ. Finally, the tips of each tree are labeled with the percentage
#' of cells assigned to each subclone.
#'
#' @examples
#' # Simulate read counts
#' sims.out<-simulate_data(N=15, S=2, M=5, alpha=0.1, beta=1.0, kappa=1, tau=999,
#'                        Ktrue=4, b.mindepth=30, b.maxdepth=50, sc.mindepth=80,
#'                        sc.maxdepth=120, scale=300, seed=8675309)
#'
#' N <- ncol(sims.out$Rs) # Number of cells
#' S <- ncol(sims.out$Rb) # Number of bulk samples
#' M <- nrow(sims.out$Rs) # Number of point mutations
#'
#' # Generate a tree
#' K<-4
#' text = paste0(paste0(paste(paste0("(", 1:(K - 1), ","),
#'              collapse = ""), K), paste(rep(")", (K - 1)), collapse = ""), ";")
#' tree <- ape::read.tree(text = text); rm(text)
#'
#' # Generate point mutations (single-nucleotide varaints) along the tree branches
#' tree$snv<-initialsnv(tree, rownames(sims.out$Rs))
#'
#' # Get the Z matrix from tree and snv
#' tree$Z <- getZ(tree)
#'
#' # Generate the proportion matrices
#' # Bulk samples
#' Pb<-t(DirichletReg::rdirichlet(S, alpha=rep(1/K,K)))
#' rownames(Pb)<-paste0('clone',1:K)
#' colnames(Pb)<-colnames(sims.out$Rb)
#'
#' # Single cells
#' Ps<-stats::rmultinom(N, 1, prob = rep(1/K, K))
#' rownames(Ps)<-paste0('clone',1:K)
#' colnames(Ps)<-colnames(sims.out$Rs)
#' tree$Pb<-Pb; rm(Pb)
#' tree$Ps<-Ps; rm(Ps)
#'
#' # Plot Z
#' plot_Z(tree)
#'
#' @export

plot_Z <- function(tree) {

  if (!inherits(tree, "phylo")){
    stop("tree must be of class \"phylo\"")
  }
  if (is.null(tree$snv)){
    stop("tree must contains point mutations (SNVs) under slot tree$snv")
  }
  if (is.null(tree$Z)){
    stop("tree must contain clonal configuration matrix, Z, under slot tree$Z")
  }
  if (is.null(tree$Ps)){
    stop("tree must contain cell-to-clone assignment matrix, Ps, under slot tree$Ps")
  }

  node_total <- max(tree$edge)
  node_shown <- ncol(tree$Z)
  node_hidden <- node_total - node_shown
  rowsums_Ps <- rowSums(tree$Ps)
  tree$tip.label[seq_len(node_shown)] <- paste0(
    "C", seq_len(node_shown),
    ": ", round(rowsums_Ps/sum(rowsums_Ps) * 100, digits = 0), "%")
  mut_id_all <- tree$Z %*% (2**seq(ncol(tree$Z), 1))
  mut_id_all <- seq(length(unique(mut_id_all)), 1)[as.factor(mut_id_all)]

  count_mut_ids<-0
  branch_ids <- rep(NA,node_total)

  # alternate binding of vectors
  trunk_ids<-(node_shown+1):node_total
  stem_ids<-1:node_shown
  alt_inds<-c(trunk_ids, stem_ids)[order(c(seq_along(trunk_ids)*2 - 1,
                                           seq_along(stem_ids)*2))]

  for (i in alt_inds) {
    mut_num <- sum(tree$snv[, 3] == i)
    if (mut_num == 0) {
      if (i == node_shown + 1) {
        branch_ids[i] <- "Root"
      } else {
        branch_ids[i] <- ""
      }
    } else{
      count_mut_ids<-count_mut_ids+1
      branch_ids[i] <-
        paste0("M", count_mut_ids, ": ", mut_num, " SNVs")
    }
  }

  pt <- ggtree::ggtree(tree)
  pt <- pt + ggplot2::geom_label(ggplot2::aes_string(x = "branch"),
                                 label = branch_ids, color = "firebrick"
  )

  pt <- pt + ggtree::geom_tiplab(hjust = 0.39, vjust = 0.85) +
    ggplot2::scale_x_reverse() + ggplot2::coord_flip()

  #top, right, bottom, left
  pt+ggplot2::theme(plot.margin = grid::unit(c(0,2,0,0), "cm"))
}

