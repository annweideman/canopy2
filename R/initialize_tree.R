#' Internal function
#'
#' Initialize phylogenetic tree.
#'
#' @examples \dontrun{
#' initialize_tree(seedling=8675309, K=4)
#' }
#' @keywords internal
initialize_tree<-function(seedling, K, Rs, S, Rb, N, Xb, Xs, alpha, tau){

  set.seed(seedling)

  # Generate a tree
  text = paste0(paste0(paste(paste0("(", 1:(K - 1), ","),
                             collapse = ""), K), paste(rep(")", (K - 1)), collapse = ""), ";")
  tree <- ape::read.tree(text = text); rm(text)
  tree$edge # The edge of the tree

  # Generate point mutations (SNV: single-nucleotide alteration) along the tree branches
  tree$snv=initialsnv(tree, rownames(Rs))

  # Get the Z matrix from tree and snv
  tree$Z = getZ(tree)
  tree$Z # M mutations x K clones

  # Generate the proportion matrix
  # Bulk samples
  Pb=t(DirichletReg::rdirichlet(S, alpha=rep(1/K,K)))
  rownames(Pb)=paste0('clone',1:K)
  colnames(Pb)=colnames(Rb)
  apply(Pb,2,sum)

  # Single cells
  Ps=stats::rmultinom(N, 1, prob = rep(1/K, K))
  rownames(Ps)=paste0('clone',1:K)
  colnames(Ps)=colnames(Rs)
  apply(Ps,2,sum)
  tree$Pb=Pb; rm(Pb)
  tree$Ps=Ps; rm(Ps)

  # Evaluate the posterior (with non-informative prior)
  tree$Post=getPost(tree, Rb, Xb, Rs, Xs, alpha, beta, kappa, tau)

  return(tree)
}
