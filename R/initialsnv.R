#' Generate SNVs
#'
#' Generate point mutations (single-nucleotide variants, SNVs) along the tree
#' branches.
#'
#' @param tree an object of class \code{pyhlo} representing the phylogenetic
#' tree.
#' @param snv.name an object of class \code{character} containing the names
#' of the mutations.
#'
#' @return
#' A matrix with five columns: \code{snv} (the SNV number), \code{snv.strt.node}
#' (the start node number), \code{snv.end.node} (the end node number), and
#' \code{snv.edge} (the edge number). Rownames correspond to mutations.
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
#' tree$snv
#'
#' @export

initialsnv <- function (tree, snv.name) {

  if (!inherits(tree, "phylo")){
    stop("tree must be of class \"phylo\"")
  }
  if (!is.character(snv.name)){
    stop("snv.name must be of class \"character\"")
  }

  # Number of point mutations
  snv.no <- length(snv.name)
  # If tree has more than two edges, then randomly place SNVs along these edges
  if (nrow(tree$edge) > 2) {
    snv.edge <- sample(2:nrow(tree$edge), size = snv.no, replace = T)
    # Otherwise, put SNVs on the only other edge besides the edge corresponding to
    # the normal clone
  } else {
    snv.edge <- rep(2, snv.no)
  }
  snv.mat <- cbind(1:snv.no, tree$edge[snv.edge, ], snv.edge)
  colnames(snv.mat) <- c("snv", "snv.strt.node", "snv.end.node", "snv.edge")
  rownames(snv.mat) <- snv.name
  return(snv.mat)
}
