#' Generate the clonal configuration matrix, Z
#'
#' Derive the binary clonal configuration matrix, Z, from a phylogenetic
#' tree.
#'
#' @param tree an object of class \code{pyhlo} representing the phylogenetic
#' tree.
#'
#' @return
#' A M (mutations) x K (clones) binary clonal configuration matrix, Z, that
#' indicates which mutations belong to which clones.
#'
#' @examples
#' # Generate a tree
#' K<-3
#' text = paste0(paste0(paste(paste0("(", 1:(K - 1), ","),
#'                            collapse = ""), K), paste(rep(")", (K - 1)), collapse = ""), ";")
#' tree <- ape::read.tree(text = text); rm(text)
#'
#' # Generate point mutations (single-nucleotide varaints) along the tree branches
#' tree$snv<-initialsnv(tree, snv.name=paste("mut",1:5))
#'
#' # Get the Z matrix from tree and SNVs
#' tree$Z <- getZ(tree)
#'
#' @export

getZ <- function(tree) {

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

  K <- (nrow(tree$edge) + 2)/2 #Number of clones
  M <- nrow(tree$snv) #Number of mutations

  # Create a matrix of zeros
  Z <- matrix(nrow = M, ncol = K, data = 0)
  rownames(Z) <- rownames(tree$snv)
  colnames(Z) <- paste("clone", 1:K, sep = "")

  clonal.snv <- vector("list", K)
  for (tip in 2:K) {

    # Set child node to be the tip (clone) of interest
    child.node <- tip
    # Find parent node corresponding to child node
    parent.node <- tree$edge[which(tree$edge[, 2] == child.node), 1]

    # As long as the parent node has value greater than the number of clones
    while (parent.node > K) {

      # Find SNVs where parent and child node correspond to start and end node
      snvtemp <- intersect(which(tree$snv[, 2] == parent.node),
                           which(tree$snv[, 3] == child.node))

      # If an SNV exists along that branch, store in appropriate list
      if (length(snvtemp) > 0) {
        clonal.snv[[tip]] <- c(clonal.snv[[tip]], snvtemp)
      }

      # Set child node equal to parent node
      child.node <- parent.node
      # If child node has reached the top of the tree, then end while loop
      if (child.node == (K + 1))
        break

      # Find parent node corresponding to updated child node
      parent.node <- tree$edge[which(tree$edge[, 2] == child.node), 1]
    }
  }
  clonal.snv[[1]] <- 0
  for (ki in 2:K) {
    # Impute value of 1 for positions where SNVs exist
    Z[clonal.snv[[ki]],ki] <- 1
  }
  return(Z)

}
