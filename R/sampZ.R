#' Internal function
#'
#' Sample the clonal configuration matrix, Z, using Metropolis Hastings (an
#' accept reject algorithm). See Alg. 1 in the main text for further details.
#'
#' @examples \dontrun{
#' sim.tree<-initialize_tree(seedling=8675309)
#' sampZ(tree=sim.tree)
#' }
#' @keywords internal

sampZ=function(tree, Rb, Xb, Rs, Xs, alpha, beta, kappa, tau, ...){
  accept<-c()
  M = nrow(tree$snv) # number of mutations
  for(i in 1:M){ # Update each mutation through a loop
    snv.new = tree$snv
    # sample from edges excluding the leftmost and the current edge
    # the probability for each edge to be selected is the same (i.e., symmetric)
    snv.edge.i = sample(setdiff(2:nrow(tree$edge), tree$snv[i,4]), size = 1)
    snv.new[i, 2:4] = c(tree$edge[snv.edge.i, 1:2], snv.edge.i)
    tree.new=tree
    tree.new$snv=snv.new
    tree.new$Z=getZ(tree.new)

    tree.new$Post=getPost(tree.new, Rb, Xb, Rs, Xs, alpha, beta, kappa, tau)
    r=exp(tree.new$Post-tree$Post)
    if(r>= stats::runif(1)){
      tree=tree.new
      accept<-c(accept,1)
    } else{
      tree=tree
      accept<-c(accept,0)
    }
  }
  return(list(tree,accept))
}
