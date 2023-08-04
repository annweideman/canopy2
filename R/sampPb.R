#' Internal function
#'
#' Sample the bulk sample to clone assignment matrix, Pb, using Metropolis
#' Hastings (an accept reject algorithm). See Alg. 1 in the main text for
#' further details.
#'
#' @examples \dontrun{
#' sim.tree<-initialize_tree(seedling=8675309)
#' sampPb(tree=sim.tree)
#' }
#' @keywords internal

sampPb=function(tree, Rb, Xb, Rs, Xs, alpha, beta, kappa, tau, ...){
  accept<-c()
  S=ncol(tree$Pb) # number of bulk samples
  K=nrow(tree$Pb) # number of clones
  for(t in 1:S){ # Update each bulk sample through a loop
    for(k in 1:K){
      Pb.new=tree$Pb
      # uniform distribution centered at the current proportion
      Pb.new[k,t]=stats::runif(1, min = max(0,Pb.new[k,t]-0.1), max=min(1,Pb.new[k,t]+0.1))
      Pb.new[,t]=Pb.new[,t]/sum(Pb.new[,t])
      tree.new=tree
      tree.new$Pb=Pb.new
      tree.new$Post=getPost(tree.new, Rb, Xb, Rs, Xs, alpha, beta, kappa, tau)
      r=exp(tree.new$Post-tree$Post)
      if(r>= stats::runif(1)){
        tree=tree.new
        accept<-c(accept,1)
      }else{
        tree=tree
        accept<-c(accept,0)
      }
    }
  }
  return(list(tree,accept))
}
