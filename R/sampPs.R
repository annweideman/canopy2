#' Internal function
#'
#' Sample the single-cell to clone assignment matrix, Ps, using Metropolis
#' Hastings (an accept reject algorithm). See Alg. 1 in the main text for
#' further details.
#'
#' @keywords internal

sampPs=function(tree, Rb, Xb, Rs, Xs, alpha, beta, kappa, tau, ...){
  accept<-c()
  N=ncol(tree$Ps) # number of cells
  K=nrow(tree$Ps) # number of clones
  for(j in 1:N){ # Update each cell through a loop
    Ps.new=tree$Ps
    # sample from clones excluding the current clone
    # the probability for each clone to be selected is the same (i.e., symmetric)
    Ps.new[,j]=stats::rmultinom(1,1,prob=(1-tree$Ps[,j])/(K-1))
    tree.new=tree
    tree.new$Ps=Ps.new
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
