#' Diagnostic Plots
#'
#' Plot posterior densities, trace, and lag autocorrelation factor (lag ACF) 
#' associated with each number of subclones. 
#' Requires output from function \code{get_trees()}.
#'
#' @param get.trees.out sample of phylogenetic trees output from function
#' \code{get_trees()} of class \code{get_trees}.
#' @param project a string specifying the project name to append to the filenames 
#' when saving output. The final filenames for \code{get_diagnostics} output will 
#' be \"project_lagACF.pdf\" (lag ACF plots), \"project_posteriors.pdf\" 
#' (posterior densities), and \"project_trace.pdf\" (trace plots). 
#' @param outpath a string specifying the location at which to save output 
#' generated from \code{get_diagnostics}.
#'
#' @return
#' @export
#'
#' @examples
#'
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

  if(dev.cur() > 1) grDevices::dev.off()
  
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
    plot.title <- ggplot2::ggtitle("Joint posterior distributions",
                                   "(medians and 95% HPD intervals)")
    chain.mat<-matrix(unlist(posteriors),ncol=nchains)
    colnames(chain.mat)<-chain.names
    bayesplot::mcmc_areas(chain.mat, pars=chain.names, prob = 0.95,
               area_method="scaled height") +
               ggplot2::labs(title = paste0("Joint posterior densities for K=",K),
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
      do.call(gridExtra::grid.arrange, c(lapply(i:length(p1), function(x) p1[[x]]),
                                list(nrow=2,ncol=2)))
    }
    else{
      do.call(gridExtra::grid.arrange, list(p1[[i]],p1[[i+1]],p1[[i+2]],p1[[i+3]],
                                   nrow=2,ncol=2))
    }
  }
  if(is.null(outpath)==F){grDevices::dev.off()}
  
  #-----------------------------------------------------------------------------
  # Plot trace
  #-----------------------------------------------------------------------------

  if(is.null(outpath)==F){
    grDevices::pdf(file = paste0(outpath,"/", project, "_trace.pdf"), onefile = TRUE, height=6)
  }

  counter<-1
  par(mfrow = c(2, 2), mar=c(4,4,4,4))
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
  
  counter<-1
  stop<-0
  for (K in Klist) {
    while(stop<nchains*counter){
      start<-1+stop
      stop<-ifelse(start+3<counter*nchains, start+3, counter*nchains)
      graphics::par(mfrow = c(1, stop-start+1),oma = c(2, 0, 4, 0))
      for(i in (start-nchains*(counter-1)):(stop-nchains*(counter-1))){
        stats::acf(post.list[[i]],
                   main=paste("Chain", eval(parse(text="i"))))
        graphics::mtext(side=3, line=1.5, cex=1.2,
                      bquote("Lag ACF after"~.(pburn*100)*"% burn-in"),
                      outer=T)
      graphics::mtext(side=3, line=0.5, cex=1.1, 
                      paste0("(K = ",K,")"), col="blue", outer=T)
      }
    }
    counter<-counter+1
  }
  
  if(is.null(outpath)==F){grDevices::dev.off()}
}
