#' Best Phylogenetic Tree
#'
#' Produce the optimal configuration of the phylogenetic tree by using a modified
#' version of the Bayesian information criterion (BIC).
#'
#' @param get.trees.out sample of phylogenetic trees output from function
#' \code{get_trees()} of class \code{get_trees}.
#' @param save.muts a logical indicating whether a text file containing all
#' mutations along the branches should be saved. If \code{TRUE}, must specify
#' \code{project} and \code{outpath}. Defaults to \code{FALSE}.
#' @param save.clones a logical indicating whether a text file containing the
#' mutations assigned to clones should be saved. If \code{TRUE}, must specify
#' \code{project} and \code{outpath}. Defaults to \code{FALSE}.
#' @param save.plot a logical indicating whether a plot of the best tree should be
#' saved. If \code{TRUE}, must specify \code{project} and \code{outpath}.
#' Defaults to \code{FALSE}.
#' @param project a string specifying the project name to add to the filename
#' for output generated from \code{save.muts} and \code{save.plot}. The final
#' filenames for \code{save.muts} and \code{save.plot} will be "projectname_muts.txt"
#' and "projectname_plot_tree.jpg," respectively.
#' @param outpath a string specifying the location at which to save output
#' generated from \code{save.muts} and \code{save.plot}.
#'
#' @details{
#' The Bayesian information criterion (BIC) is computed by modifying the
#' classical definition \insertCite{Schwarz_1978}{canopy2}. Since the prior is flat, the
#' likelihood is proportional to the posterior, so the maximum a posteriori (MAP)
#' estimate is equivalent to the maximum likelihood estimate (MLE). The MAP
#' estimate is equal to the mode of the posterior distribution and defined as
#' the value of the parameter that maximizes the posterior distribution, which
#' is the value at which the distribution reaches its highest peak. Thus, we can
#' replace the maximized log-likelihood with the maximized log-posterior.
#' However, since we have multiple chains, we take the mean of this maximized
#' value across \code{nchains}.
#'
#' \deqn{\text{BIC}_{\text{MAP}} = -2\log\left(\frac{1}{n_{chains}}\sum_{j=1}^{n_{chains}}p_j(\widehat{\theta}_{\text{MAP}}|y)\right)+p\log(n),}
#'
#' where \eqn{n_{chains}} denotes the number of MCMC chains,
#' \eqn{\widehat{\theta}_{\text{MAP}}} denotes the MAP estimate,
#' \eqn{y} denotes the observed data, \eqn{p} denotes the number of parameters,
#' and \eqn{n} denotes the total sample size.
#'
#' In this particular case, \eqn{p=K} (the number of subclones) since we are
#' attempting to determine the optimal number of subclones. The sample size,
#' \eqn{n = 2(M \times N) + 2(M \times T)}, reflects the aggregate dimensions of
#' the observed data where \eqn{\text{dim}(R^{s})=\text{dim}(X^{s})=M\times N} and
#' \eqn{\text{dim}(R^{b})=\text{dim}(X^{b})=M\times T} for
#' \eqn{M} mutations, \eqn{N} single cells, and \eqn{T} bulk samples.
#' }
#'
#' @return
#' A list containing: \code{K}, a numeric corresponding to the optimal number of
#' subclones occurring at the minimum BIC, \code{tree}, an object of class
#' \code{"phylo"} representing the best phylogenetic tree corresponding to the
#' optimal number of subclones \code{K}, \code{branches.muts}, a matrix containing the
#' names of the mutations assigned to each branch along the tree,
#' \code{clonals.muts}, a matrix containing the names of the mutations assigned
#' to each subclone, \code{posteriors}, a numeric of posteriors corresponding to
#' the final tree, \code{BIC}, a numeric corresponding to the minimum BIC across
#' all subclones, and \code{acceptance.rate}, a numeric corresponding to the MCMC
#' acceptance rate for the best tree.
#' @import magrittr
#'
#' @examples
#' #Load post-processed data for patient GBM10
#' data("GBM10_postproc")
#'
#' # Run Canopy2 to get list of phylogenetic trees corresponding to all chains
#' # and all subclones.
#' # Note: this is written to quickly compile for demonstration purposes, as we
#' # know the optimal number of subclones is 6. In practice, we would attempt a
#' # larger range of subclones (e.g., Klist=3:10) and a larger number of
#' # iterations and chains (e.g., niter=50000 and nchains=10).
#' get.trees.out<-get_trees(Rs=GBM10_postproc@Rs, Rb=GBM10_postproc@Rb,
#'                          Xs=GBM10_postproc@Xs, Xb=GBM10_postproc@Xb,
#'                          alpha=GBM10_postproc@param.est$alpha,
#'                          beta=GBM10_postproc@param.est$beta, kappa=1,
#'                          tau=999, Klist=5:7, niter=5000, nchains=5, thin=20,
#'                          pburn=0.2, seed=8675309)
#'
#' # Examine diagnostic plots
#' get_diagnostics(get.trees.out)
#'
#' # Get best tree across all chains and subclones via BIC
#' best.tree.out<-get_best_tree(get.trees.out)
#'
#' best.tree.out
#'
#' @export

get_best_tree<-function(get.trees.out,
                        save.muts=F,
                        save.clones=F,
                        save.plot=F,
                        project=NULL,
                        outpath=NULL){

  # check arguments
  if (!inherits(get.trees.out, "get_trees")){
    stop("get.trees.out must be of class \"get_trees\" - output from the get_trees
         function")
  }
  if(!is.logical(save.muts)){
    stop("Argument save.muts must be a logical (TRUE or FALSE).")
  }
  if(!is.logical(save.clones)){
    stop("Argument save.clones must be a logical (TRUE or FALSE).")
  }
  if(!is.logical(save.plot)){
    stop("Argument save.plot must be a logical (TRUE or FALSE).")
  }
  if(is.null(project) & (save.muts | save.clones | save.plot)){
    stop("You must specify a project name in order to save the output.")
  }
  if(!is.null(project) & !is.character(project)){
    stop("Argument project must be of type string.")
  }
  if(!is.null(outpath) & is.null(project)){
    stop("Argument project must be specified if argument outpath is specified.")
  }
  if(is.null(outpath) & !is.null(project)){
    stop("Argument outpath must be specified if argument project is specified.")
  }
  if(!is.null(outpath) & !(save.muts | save.clones | save.plot)){
    stop("You must set save.muts, save.clones, and/or save.plot to 'TRUE' so
         that the output can be directed to the location specified by \"outpath\".")
  }
  if(is.null(outpath) & (save.muts | save.clones | save.plot)){
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
  Klist<-get.trees.out$Klist # Range of subclones to evaluate
  nchains<-get.trees.out$nchains # Number of chains of MCMC

  # Initialize
  final.out<-list()
  counter<-0

  # Loop over subclones in pre-specified range Klist
  for (K in Klist){

    # Apply kernel density estimation (KDE) to the log-posterior estimates, and
    # then find the maximum of this kernel density. This would represent the
    # log-posterior evaluated at the maximum a posteriori (MAP) estimate.
    density_at_mode <- function(data) {

       # Perform a kernel density estimate
       densities <- stats::density(data)

       # Return value of log-posterior at this maximum density
       densities$x[which.max(densities$y)]
    }

    # Evaluate the log-posterior at its MAP for each chain
    list.maxLogPost<-lapply(1:nchains, function(x)
                     density_at_mode(samples[[x+counter]]$posteriors))

    # Store ID associated with the best chain which maximizes the log-posterior
    id.max.post<-which.max(list.maxLogPost)

    # Store samples associated with the best chain
    best.chain<-samples[[id.max.post+counter]]

    # Log-posteriors from best chain
    logPost<-best.chain$posteriors

    # Compute the mean of the maximized log-Posteriors across all chains
    mean.maxLogPost<-mean(unlist(list.maxLogPost))

    # Compute BIC
    BIC<- -2*mean.maxLogPost + K*log(2*M*N+2*M*S)

    # Store the final tree from the best chain
    final.tree<-best.chain$tree[[length(best.chain$tree)]]

    # Store output
    temp.out<-list("K"=K,
                   "tree"=final.tree,
                   "Ps"=final.tree$Ps,
                   "posteriors"=best.chain$posteriors,
                   "BIC"=BIC,
                   "acceptance.rate"=mean(best.chain$accept),
                   "mean.maxLogpost"=mean.maxLogPost)
    final.out<-append(final.out,temp.out)

    # Update counter to skip to next K (each K has nchains)
    counter<-counter+nchains
  }

  # Generate lists of each parameter of interest
  Klist<-sapply(seq(1,length(final.out),7), function(x) final.out[x])
  tree.list<-sapply(seq(2,length(final.out),7), function(x) final.out[x])
  Ps.list<-sapply(seq(3,length(final.out),7), function(x) final.out[x])
  post.list<-sapply(seq(4,length(final.out),7), function(x) final.out[x])
  BIC.list<-sapply(seq(5,length(final.out),7), function(x) final.out[x])
  accept.list<-sapply(seq(6,length(final.out),7), function(x) final.out[x])
  mean.maxLogPost.list<-sapply(seq(7,length(final.out),7), function(x) final.out[x])

  # Print BIC update to console
  for(i in 1:length(Klist)){
    print(paste0("k=", Klist[[i]], "; mean of maximized log-posteriors across all chains=",
                 round(mean.maxLogPost.list[[i]],2),
                 "; BIC=", round(BIC.list[i]$BIC,2)))
  }

  # Locate the list ID for which the minimum BIC occurs
  BIC.list<-unlist(BIC.list)
  id.min.BIC<-which.min(BIC.list)

  # Store posteriors at minimum BIC
  post.min.BIC<-post.list[id.min.BIC]$posteriors

  # Store K (number of subclones) at minimum BIC
  K.min.BIC<-Klist[id.min.BIC]$K

  # Store tree at minimum BIC
  tree.min.BIC<-tree.list[id.min.BIC]$tree

  # Store minimum BIC
  min.BIC<-BIC.list[id.min.BIC]

  # Plot BIC vs number of subclones to visually locate optimal K (occurs at min BIC)
  plot(Klist, BIC.list,
       xlab = "Number of subclones",
       ylab = "BIC",
       type = "b")
  graphics::axis(1, at = K)
  graphics::abline(v = Klist[which.min(BIC.list)], lty = 2)
  graphics::title(paste("Optimal number of subclones for Canopy2 (min BIC)"))

  # Acceptance rate at minimum BIC
  accept.min.BIC<-accept.list[id.min.BIC]$acceptance.rate

  # Plot best inferred phylogeny and output mutations at each node
  branches.muts<-plot_tree(tree=tree.min.BIC,
                           save.muts=save.muts,
                           save.plot=save.plot,
                           outpath=outpath,
                           project=project)

  # Return the single nucleotide variants (SNVs) belonging to each subclone
  clonal.muts<-get_clonal_composition(tree=tree.min.BIC)

  # if TRUE, save text file containing all mutations along the branches
  if (save.muts){
    utils::write.table(branches.muts, file = paste0(outpath,"/", project, "_branch_muts.txt"),
                       col.names = FALSE, row.names = FALSE,
                       quote = FALSE, sep = '\t')
  }

  # if TRUE, save text file containing mutations belonging to each clone
  if (save.clones){
    utils::write.table(clonal.muts, file = paste0(outpath,"/", project, "_clonal_muts.txt"),
                       col.names = FALSE, row.names = FALSE,
                       quote = FALSE, sep = '\t')
  }

  return(list("K"=K.min.BIC, # Subclones at minimum BIC
              "tree"=tree.min.BIC, # Tree at minimum BIC
              "branches.muts"=branches.muts, # mutations within each cluster along tree branches
              "clonal.muts"=clonal.muts, # mutations belonging to each subclone
              "posteriors"=post.min.BIC, # Posteriors for final tree
              "BIC"=min.BIC, # Minimum BIC
              "acceptance.rate"=accept.min.BIC #MCMC acceptance rate
              )
         )
}
