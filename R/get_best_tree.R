#' Best Phylogenetic Tree
#'
#' Produce the optimal configuration of the phylogenetic tree by using the 
#' deviance information criterion (DIC). 
#'
#' @param get.trees.out sample of phylogenetic trees output from function
#' \code{get_trees()} of class \code{get_trees}.
#' @param method type of DIC used for model selection and must be
#' one of method \code{"spiegelhalter"} or \code{"gelman"}. Defaults to
#' \code{"spiegelhalter"}.
#' @param save.muts a logical indicating whether a text file containing all
#' mutations along the branches should be saved. If \code{TRUE}, must specify
#' \code{project} and \code{outpath}. Defaults to \code{FALSE}.
#' @param save.plot a logical indicating whether plot of the best tree should be 
#' saved. If \code{TRUE}, must specify \code{project} and \code{outpath}. 
#' Defaults to \code{FALSE}.
#' @param project a string specifying the project name to add to the filename 
#' for output generated from \code{save.muts} and \code{save.plot}. The final 
#' filenames for \code{save.muts} and \code{save.plot} will be project_muts.txt 
#' and project_plot_tree.jpg, respectively.
#' @param outpath a string specifying the location at which to save output 
#' generated from \code{save.muts} and \code{save.plot}.
#'
#' @return
#' @export
#'
#' @examples
get_best_tree<-function(get.trees.out,
                        method="spiegelhalter", 
                        save.muts=F,
                        save.plot=F,
                        project=NULL, 
                        outpath=NULL){
  
  # check arguments
  if (!inherits(get.trees.out, "get_trees")){
    stop("get.trees.out must be of class \"get_trees\" - output from the get_trees
         function")
  }
  if (!(method %in% c("spiegelhalter","gelman"))){
    stop("method must be either 'spiegelhalter' or 'gelman'")
  }
  if(is.null(project) & (save.muts | save.plot)){
    stop("You must specify a project name in order to save the output.")
  }
  if(!is.null(project) & !is.character(project)){
    stop("The project argument must be of type string.")
  }
  if(!is.null(outpath) & is.null(project)){
    stop("Argument project must be specified if argument outpath is specified.")
  }
  if(is.null(outpath) & !is.null(project)){
    stop("Argument outpath must be specified if argument project is specified.")
  }
  if(!is.null(outpath) & !(save.muts | save.plot)){
    stop("You must set save.muts and/or save.plot to 'TRUE' so that the output
         can be directed to the location specified by \"outpath\".")
  }
  if(is.null(outpath) & (save.muts | save.plot)){
    stop("You must specify an outpath in order to save the output.")
  }
  if(!is.null(outpath) & !is.character(outpath)){
    stop("The outpath argument must be of type string.")
  }
  if(!is.null(outpath)){
    if(!file.exists(outpath)){
      stop("The location specified for outpath is not valid.")
    }
  }
  
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

  # Initialize
  final.out<-list()
  counter<-0
  
  # Loop over subclones in pre-specified range, Klist
  for (K in Klist){
    
    # Generate list of mean log-posteriors for each chain
    list.mean.post<-lapply(1:nchains, function(x) 
                           mean(samples[[x+counter]]$posteriors))
    
    # Store ID associated with the best chain which maximizes the log-posterior
    id.max.post<-which.max(list.mean.post)
    
    # Store samples associated with the best chain
    best.chain<-samples[[id.max.post+counter]]
    
    # Compute the mean of the posterior samples from the best chain
    logPost.mean<-mean(best.chain$posteriors)

    # Compute DIC under Gelman definition
    if(method=="gelman"){
      
      # Compute DIC via Gelman defn
      DIC<-.5*var(best.chain$posteriors) + -2*logPost.mean
      
    } # End if
    
    # Compute DIC via Spiegelhalter defn    
    if(method=="spiegelhalter"){

      # Parameter means (theta_bar) for DIC
      Z.list<-lapply(1:length(best.chain$tree),
                     function(x) best.chain$tree[[x]]$Z)
      Z.mean<-apply(simplify2array(Z.list), 1:2, mean)

      Ps.list<-lapply(1:length(best.chain$tree),
                      function(x) best.chain$tree[[x]]$Ps)
      Ps.mean<-apply(simplify2array(Ps.list), 1:2, mean)

      Pb.list<-lapply(1:length(best.chain$tree),
                      function(x) best.chain$tree[[x]]$Pb)
      Pb.mean<-apply(simplify2array(Pb.list), 1:2, mean)

      # Compute DIC 
      DIC<-dic_fun(logPost.mean, Z.mean, Ps.mean,
                   Pb.mean, Rb, Xb, Rs, Xs, alpha, beta, kappa, tau)
    } # End if  

    # Store the id of the final tree corresponding to best chain
    id.tree<-length(best.chain$posteriors)

    # Store the final tree
    final.tree<-best.chain$tree[[id.tree]]

    # Sort Ps to aggregate 1's along main diagonal for better visualization
    temp.mat<-matrix(NA,nrow=nrow(final.tree$Ps),ncol=1)
    rownames(temp.mat)<-rownames(final.tree$Ps)
    id.ones<-lapply(1:nrow(final.tree$Ps),function(x) which(final.tree$Ps[x,]==1))

    # Generate mini binary matrices that are concatenated to create final matrix
    for(list.id in length(id.ones):1){
        temp.mat0<-matrix(0,nrow=nrow(final.tree$Ps),ncol=length(id.ones[[list.id]]))
        temp.mat0[list.id,]<-1
        colnames(temp.mat0)<-names(id.ones[[list.id]])
        temp.mat<-cbind(temp.mat,temp.mat0)
    }
     
    # Remove vector of NAs that was introduced in previous step
    temp.mat<-temp.mat[,-1]
    # Replace Ps in final tree with new Ps that has 1's along main diagonal
    final.tree$Ps<-temp.mat
  
    # Store output
    temp.out<-list("K"=K,
                   "tree"=final.tree,
                   "Ps"=final.tree$Ps,
                   "posteriors"=best.chain$posteriors,
                   "DIC"=DIC,
                   "acceptance.rate"=mean(best.chain$accept))
    final.out<-append(final.out,temp.out)
    
    # Update counter to skip to next K (each K has nchains)
    counter<-counter+nchains
  }
  
  # Generate lists of each parameter of interest
  Klist<-sapply(seq(1,length(final.out),6), function(x) final.out[x])
  tree.list<-sapply(seq(2,length(final.out),6), function(x) final.out[x])
  Ps.list<-sapply(seq(3,length(final.out),6), function(x) final.out[x])
  post.list<-sapply(seq(4,length(final.out),6), function(x) final.out[x])
  DIC.list<-sapply(seq(5,length(final.out),6), function(x) final.out[x])
  accept.list<-sapply(seq(6,length(final.out),6), function(x) final.out[x])
  
  # Print DIC update to console
  for(i in 1:length(Klist)){
    print(paste0("k=", Klist[i], "; posterior mean from best chain=",
                 round(mean(post.list$posteriors[i]),2), 
                 "; DIC=", round(DIC.list[i]$DIC,2)))
  } # End for loop
  
  # Locate the list ID for which the minimum DIC occurs
  DIC.list<-unlist(DIC.list)
  id.min.DIC<-which.min(DIC.list)
  
  # Store K (number of subclones) at minimum DIC
  K.min.DIC<-Klist[id.min.DIC]$K
  
  # Store best tree at minimum DIC
  tree.min.DIC<-tree.list[id.min.DIC]$tree
  
  # Store posteriors at minimum DIC
  post.min.DIC<-post.list[id.min.DIC]$posteriors
  
  # Store minimum DIC
  min.DIC<-DIC.list[id.min.DIC]
  
  graphics.off()
  
  # Plot DIC vs number of subclones to visually locate optimal K (occurs at min DIC)
  plot(Klist, DIC.list, 
       xlab = "Number of subclones", 
       ylab = "DIC", 
       type = "b", 
       xaxt = "n")
  axis(1, at = K)
  abline(v = Klist[which.min(DIC.list)], lty = 2)
  title(paste("Optimal number of subclones for Canopy2 (min DIC)"))
  
  # Acceptance rate at minimum DIC
  accept.min.DIC<-accept.list[id.min.DIC]$acceptance.rate
  
  # Plot best inferred phylogeny and output mutations at each node
  inferred.mut.list<-plot_tree(tree=tree.min.DIC,
                               save.muts=save.muts,
                               save.plot=save.plot, 
                               outpath=outpath, 
                               project=project)

  # Return the single nucleotide variants (SNVs) belonging to each subclone
  clonal.mut.list<-get_clonal_composition(tree=tree.min.DIC)
  names(clonal.mut.list)<-paste0("clone",tree.min.DIC$tip.label)
  
  return(list("K"=K.min.DIC, # Subclones at minimum DIC
              "tree"=tree.min.DIC, # Tree at minimum DIC
              "mut.list"=inferred.mut.list, # List of mutations within each cluster along tree branches
              "clonal.mut.list"=clonal.mut.list, # List of mutations belonging to each subclone
              "posteriors"=post.min.DIC, # Posteriors for final tree
              "DIC"=min.DIC, # Minimum DIC
              "acceptance.rate"=accept.min.DIC #MCMC acceptance rate
              ) 
         )
}
