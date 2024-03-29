#' Diagnostic Plots
#'
#' Plot posterior densities, trace, and lag autocorrelation factor (lag ACF)
#' associated with each number of subclones, K.
#' Requires output from function \code{get_trees()}.
#'
#' @param get.trees.out sample of phylogenetic trees output from function
#' \code{get_trees()} of class \code{get_trees}.
#' @param project a string specifying the project name to append to the filenames
#' when saving output. The final filenames for \code{get_diagnostics} output will
#' be "projectname_lagACF.pdf" (lag ACF plots), "projectname_posteriors.pdf"
#' (posterior densities), and "projectname_trace.pdf" (trace plots).
#' @param outpath a string specifying the location at which to save output
#' generated from \code{get_diagnostics}.
#'
#' @details
#' Diagnostic plots are helpful for determining if there is evidence of
#' convergence within the MCMC chains and also to assess agreement of these
#' plots with the tree selected as opitmal from \code{\link[canopy2]{get_best_tree}}.
#'
#' The posterior densities from the individual chains should ideally sit close
#' in space and each of the chains in the trace plots should be indistinguishable
#' except for random noise. The lag autocorrelation factor (lag ACF) should drop
#' to zero rather rapidly; large spikes that persist after lag 5 or 10 can be
#' indicative of strong within chain correlation (autocorrelation). A large
#' autocorrelation can imply that the sampling space is not properly explored
#' and can result in convergence to a local (rather than global) optimum or
#' convergence failure, altogether.
#'
#' @return
#' For each K, generates graphical displays of the posterior densities, trace,
#' and lag ACF. These graphs can be viewed in R or exported to pdfs.
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
#' get_diagnostics(get.trees.out=get.trees.out)
#'
#' @export

get_diagnostics<-function(get.trees.out, project=NULL, outpath=NULL){

  if (!inherits(get.trees.out, "get_trees")){
    stop("get.trees.out must be of class \"get_trees\" output from the get_trees
         function")
  }
  if(!is.null(project) & !is.character(project)){
    stop("The project argument must be of type string.")
  }
  if(!is.null(outpath) & !is.character(outpath)){
    stop("The outpath argument must be of type string.")
  }
  if(!is.null(outpath)){
    if(!file.exists(outpath)){
      stop("The location specified for outpath is not valid.")
    }
  }

  #if(grDevices::dev.cur() > 1) grDevices::dev.off()

  samples<-get.trees.out$samples
  nchains<-get.trees.out$nchains
  Klist<-get.trees.out$Klist
  pburn<-get.trees.out$pburn

  #-----------------------------------------------------------------------------
  # Plot posterior densities
  #-----------------------------------------------------------------------------
  post.list<-lapply(1:length(samples), function(x) samples[[x]]$posteriors)

  # Function to plot posterior densities
  plot.posteriors<-function(posteriors, K, nchains){

    chain.names<-sapply(1:nchains, function(i) paste("Chain",i))
    chain.mat<-matrix(unlist(posteriors),ncol=nchains)
    colnames(chain.mat)<-chain.names
    bayesplot::mcmc_areas(chain.mat, pars=chain.names, prob = 0.95,
                          area_method="scaled height") +
      ggplot2::labs(title = paste0("Joint log-posterior densities for K=",K),
                    subtitle = "(medians and 95% HPD intervals)")+
      ggplot2::theme(plot.title = ggplot2::element_text(size=13,hjust=0.5),
                     plot.subtitle = ggplot2::element_text(hjust=0.5),
                     plot.caption = ggplot2::element_text(color="blue", hjust=0.5))

  }

  # Arrange output in a 2x2 grid
  p1<-list()
  i<-1
  counter<-1
  for(K in Klist){
    p1[[i]]<-plot.posteriors(post.list[counter:(counter+nchains-1)], K, nchains)
    i<-i+1
    counter<-counter+nchains
  }

  if(!is.null(outpath)){
    grDevices::pdf(file = paste0(outpath,"/", project, "_posteriors.pdf"),onefile = TRUE)
  }

  iters<-seq(1,length(p1),4)
  for (i in iters) {
    if(i==max(iters)){
      p2<-do.call(gridExtra::grid.arrange, c(lapply(i:length(p1), function(x) p1[[x]]),
                                             list(nrow=2,ncol=2)))
    }
    else{
      p2<-do.call(gridExtra::grid.arrange, list(p1[[i]],p1[[i+1]],p1[[i+2]],p1[[i+3]],
                                                nrow=2,ncol=2))
    }
    p2
  }
  if(is.null(outpath)==F){grDevices::dev.off()}

  #-----------------------------------------------------------------------------
  # Plot trace
  #-----------------------------------------------------------------------------

  if(is.null(outpath)==F){
    grDevices::pdf(file = paste0(outpath,"/", project, "_trace.pdf"), onefile = TRUE, height=6)
  }

  counter<-1
  graphics::par(mfrow = c(2, 2), mar=c(4,4,4,4))
  for (K in Klist) {
    coda::traceplot(Map(coda::as.mcmc,post.list[counter:(counter+nchains-1)]),ylab="Posterior")
    #main=bquote("Trace plot for posterior after"~.(pburn*100)*"% burn-in"))
    graphics::mtext(bquote("Trace plot for posterior after"~.(pburn*100)*"% burn-in")
                    ,side=3,line=-2,outer=TRUE)
    graphics::mtext(side=3, line=0.5, cex=1, paste0("K = ",K,""), col="blue")
    counter<-counter+nchains
  }
  if(is.null(outpath)==F){grDevices::dev.off()}

  #-----------------------------------------------------------------------------
  # Plot lag autocorrelation factor (lag ACF)
  #-----------------------------------------------------------------------------

  if(is.null(outpath)==F){
    grDevices::pdf(file = paste0(outpath,"/", project, "_lagACF.pdf"),
                   onefile = TRUE, height=5)
  }

  chain.id<-1
  count<-1
  start<-1
  end<-(Klist-(min(Klist)-1))*nchains
  if(nchains<=5){n.row=1}else{n.row=2}
  for (K in Klist){
    graphics::par(mfrow = c(n.row, 5),oma = c(2, 0, 6, 0))
    for(i in start:end[count]){
      stats::acf(post.list[[i]],
                 main=paste("Chain", eval(parse(text="chain.id"))))
      graphics::mtext(side=3, line=1.5, cex=1.2,
                      bquote("Lag ACF after"~.(pburn*100)*"% burn-in"),
                      outer=T)
      graphics::mtext(side=3, line=0.5, cex=1.1,
                      paste0("(K = ",K,")"), col="blue", outer=T)
      chain.id<-chain.id+1
    }
    start<-end[count]+1
    count<-count+1
    chain.id<-1
  }

  if(is.null(outpath)==F){grDevices::dev.off()}
}
