#' Best Phylogeny
#'
#' Produce the optimal phylogenetic tree from simulations based on a
#' user-supplied information criterion.
#'
#' @param get.trees.out sample of phylogenetic trees output from function
#' \code{get_trees()} of class \code{get_trees}.
#' @param criterion information criterion used for model selection and must be
#' one of \code{"AIC"}, \code{"BIC"}, and \code{"DIC"}.
#' @param save.muts a logical indicating whether a text file containing all
#' mutations along the branches should be saved. If \code{TRUE}, must specify
#' \code{project} and \code{outpath}.
#' @param save.plot a logical indicating whether a plot of the best tree should 
#' be saved. If \code{TRUE}, must specify \code{project} and
#' \code{outpath}.
#' @param project a project name to add to the filename for output generated
#' from \code{save.muts} and \code{save.plot}.
#' @param outpath a location at which to save output generated from
#' \code{save.muts} and \code{save.plot}.
#'
#' @return
#' @export
#'
#' @examples
get_best_tree<-function(get.trees.out, criterion=c("AIC","BIC","DIC"), 
                             save.muts=F, save.plot=F, project=NULL, 
                             outpath=NULL){
  
  # check arguments
  if (!inherits(get.trees.out, "get_trees")){
    stop("get.trees.out must be of class \"get_trees\" - output from the get_trees
         function")
  }
  if(is.null(project) & (save.muts | save.plot)){
    stop("You must specify a project name in order to save the plot and/or save the list of mutations.")
  }
  if(!is.null(project) & !is.character(project)){
    stop("The project argument must be of type string.")
  }
  if(is.null(outpath) & (save.muts | save.plot)){
    stop("You must specify an outpath in order to save the plot and/or save the list of mutations.")
  }
  if(!is.null(outpath) & !is.character(outpath)){
    stop("The outpath argument must be of type string.")
  }
  if(!is.null(outpath)){
    if(!file.exists(outpath)){
      stop("The location specified for outpath is not valid.")
    }
  }
  
  temp.out<-list()
  final.out<-list()
  counter<-0
  
  # Grab parameters from output generated from get_trees()
  samples<-get.trees.out$samples # Sample space of phylogenetic trees
  N<-get.trees.out$N # Number of cells
  S<-get.trees.out$S # Number of bulk samples
  M<-get.trees.out$M # Number of point mutations
  Rs<-get.trees.out$Rs # Single-cell alternative read counts
  Rb<-get.trees.out$Rb # Bulk alternative read counts
  Xs<-get.trees.out$Xs # Single-cell total read counts
  Xb<-get.trees.out$Xb # Bulk total read counts
  alpha<-get.trees.out$alpha # Gene activation rates
  beta<-get.trees.out$beta # Gene deactivation rates
  kappa<-get.trees.out$kappa; tau<-get.trees.out$tau # Sequencing error
  Klist<-get.trees.out$Klist # Possible numbers of subclones
  nchains<-get.trees.out$nchains # Number of chains of MCMC
  
  for (K in Klist){
    
    # Generate list of max posteriors for each chain
    list.max.post<-lapply(1:nchains, function(x) max(samples[[x+counter]]$post))
    
    # Find the chain corresponding to the max of the list of maxes
    chain.max.post<-which.max(list.max.post)
    
    # Find the posterior corresponding to the max of the list of maxes
    max.post<-list.max.post[[chain.max.post]]
    
    #Compute information criteria
    if(criterion=="BIC"){info.criterion=bic_fun(logPost=max.post, K, M, N, S)}
    if(criterion=="AIC"){info.criterion=aic_fun(logPost=max.post, K, M, N, S)}
    if(criterion=="DIC"){
      
      # parameter means (theta_bar) for DIC
      Z.list<-lapply(1:length(samples[[chain.max.post+counter]]$tree),
                     function(x) samples[[chain.max.post+counter]]$tree[[x]]$Z)
      Z.mean<-apply(simplify2array(Z.list), 1:2, mean)
      
      Ps.list<-lapply(1:length(samples[[chain.max.post+counter]]$tree),
                      function(x) samples[[chain.max.post+counter]]$tree[[x]]$Ps)
      Ps.mean<-apply(simplify2array(Ps.list), 1:2, mean)
      
      Pb.list<-lapply(1:length(samples[[chain.max.post+counter]]$tree),
                      function(x) samples[[chain.max.post+counter]]$tree[[x]]$Pb)
      Pb.mean<-apply(simplify2array(Pb.list), 1:2, mean)
      
      info.criterion=dic_fun(samples[[chain.max.post+counter]], Z.mean, Ps.mean,
                             Pb.mean, Rb, Xb, Rs, Xs, alpha, beta, kappa, tau)}
    
    # Find the within-chain id corresponding to final.post
    # If more than one (duplicates), then just grab the first
    id.max.post<-which(samples[[chain.max.post+counter]]$post==max.post)[1]
    
    # Store the final tree and Ps
    final.tree<-samples[[chain.max.post+counter]]$tree[id.max.post]
    final.Ps<-samples[[chain.max.post+counter]]$tree[[id.max.post]]$Ps
    
    # Compute acceptance rate for selected chain
    accept.rate<-mean(samples[[chain.max.post+counter]]$accept)
    
    # Store output
    temp.out<-list("K"=K,
                   "tree"=final.tree,
                   "Ps"=final.Ps,
                   "posteriors"=samples[[chain.max.post+counter]]$post,
                   "information.criterion"=info.criterion,
                   "acceptance.rate"=accept.rate)
    final.out<-append(final.out,temp.out)
    counter<-counter+nchains
  }
  
  # Generate lists of each parameter of interest
  Klist<-sapply(seq(1,length(final.out),6), function(x) final.out[x])
  tree.list<-sapply(seq(2,length(final.out),6), function(x) final.out[x])
  Ps.list<-sapply(seq(3,length(final.out),6), function(x) final.out[x])
  post.list<-sapply(seq(4,length(final.out),6), function(x) final.out[x])
  info.criterion.list<-sapply(seq(5,length(final.out),6), function(x) final.out[x])
  accept.list<-sapply(seq(6,length(final.out),6), function(x) final.out[x])
  
  # Find the list ID for which the minimum information criterion occurs
  info.criterion.list<-unlist(info.criterion.list)
  id.min.info.criterion<-which.min(info.criterion.list)
  
  # K (number of subclones) at minimum information criterion
  K.min.info.criterion<-Klist[id.min.info.criterion]$K
  
  # Tree at minimum information criterion
  tree.min.info.criterion<-tree.list[id.min.info.criterion]$tree[[1]]
  
  # Posteriors at minimum information criterion
  post.min.info.criterion<-post.list[id.min.info.criterion]$post
  
  # Minimum information criterion
  min.info.criterion<-info.criterion.list[id.min.info.criterion]
  
  # Acceptance rate at minimum information criterion
  accept.min.info.criterion<-accept.list[id.min.info.criterion]$acceptance.rate
  
  # Plot best inferred phylogeny and output mutations at each node
  inferred.mut.list<-plot_tree_test(tree=tree.min.info.criterion, Rs, N, S,
                                    ptitle=paste0("Inferred Tree: K=",
                                                  K.min.info.criterion), annovar=NULL, save.muts=save.muts,
                                    save.plot=save.plot, outpath=outpath, project=project)
  
  return(list("K"=K.min.info.criterion,
              "tree"=tree.min.info.criterion,
              "inferred mutations"=inferred.mut.list,
              "posteriors"=post.min.info.criterion,
              "information criterion"=min.info.criterion,
              "acceptance rate"=accept.min.info.criterion))
  
}
