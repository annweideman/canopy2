#' Sample space of phylogenetic trees
#'
#' Produce a list of phylogenetic trees in a sample space of size
#' (number of chains) x (number of subclones). Sampling is performed using
#' Metropolis-within-Gibbs. Output is of class \code{get_trees}.
#'
#' @param Rs \eqn{M} (mutations) x \eqn{N} (cells) matrix of single-cell
#' alternative read counts. Rownames must match those of argument \code{Rb}
#' @param Rb \eqn{M} (mutations) x \eqn{S} (bulk samples) matrix of bulk
#' alternative read counts. Rownames must match those of argument \code{Rs}
#' @param Xs \eqn{M} (mutations) x \eqn{N} (cells) matrix of single-cell
#' total (benign + mutated) read counts. Rownames must match those of argument
#' \code{Xb}.
#' @param Xb \eqn{M} (mutations) x \eqn{S} (bulk samples) matrix of bulk total
#' (benign + mutated) read counts. Rownames must match those of argument
#' \code{Xs}.
#' @param alpha numeric vector of positive values with length \eqn{M}
#' (number of mutations) representing the gene activation rates
#' @param beta numeric vector of positive values with length \eqn{M}
#' (number of mutations) representing the gene deactivation rates
#' @param kappa a positive value used in the computation of the
#' sequencing error defined as \eqn{\kappa/(\kappa+\tau)}
#' @param tau a positive value used in the computation of the sequencing error
#' defined as \eqn{\kappa/(\kappa+\tau)}
#' @param Klist numeric vector containing the possible numbers of subclones.
#' @param niter number of iterations of MCMC. Defaults to 10000.
#' @param nchains number of chains for MCMC. Defaults to 20.
#' @param thin a number representing the increment at which to store MCMC output.
#'  For example, if \code{thin=10} then every 10th iteration is stored. Defaults
#'  to 20.
#' @param pburn a decimal denoting the percentage of burn-in to store. Defaults
#' to 0.50 (50\\%).
#' @param ncores the number of cores to use for parallel execution. If not
#' specified, defaults to one-half the number of cores detected by the
#' parallelly package.
#' @param seed a state (positive integer) to set the random number generation.
#' Defaults to 8675309.
#'
#' @return
#' An object of class "get_trees" containing: \code{N}, the number of
#' single-cells, \code{S}, the number of bulk samples, \code{M}, the number of
#' mutations, \code{Rs}, the single-cell alternative read count matrix,
#' \code{Rb} the bulk alternative read count matrix, \code{Xs}, the single-cell
#' total (benign + mutated) read count matrix, \code{Xb}, the bulk total
#' (benign + mutated) read count matrix, \code{alpha}, a vector of
#' mutation-specific gene activation rates, \code{beta}, a vector
#' of mutation-specific gene deactivation rates, \code{kappa}, a numeric used to
#' calculate the sequencing error rate, \code{tau}, a numeric used to calculate
#' the sequencing error rate, \code{Klist}, a range of subclones to try,
#' \code{nchains}, the number of MCMC chains, and \code{pburn}, the percentage of
#' burn-in to remove.
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

get_trees<-function(Rs, Rb, Xs, Xb, alpha, beta, kappa, tau,
                    Klist, niter=10000, nchains=20, thin=20, pburn=0.5,
                    ncores=NULL, seed=8675309){

  # check arguments
  if (!inherits(Rs, "matrix")){
    stop("Rs must be of class \"matrix\"")
  }
  if(!all(Rs==round(Rs)) | !is.numeric(Rs) | !all(Rs>=0)){
    stop("Rs must contain positive integers")
  }
  if (!inherits(Rb, "matrix")){
    stop("Rb must be of class \"matrix\"")
  }
  if(!all(Rb==round(Rb)) | !is.numeric(Rb) | !all(Rb>=0)){
    stop("Rb must contain positive integers")
  }
  if (!inherits(Xs, "matrix")){
    stop("Xs must be of class \"matrix\"")
  }
  if(!all(Xs==round(Xs)) | !is.numeric(Xs) | !all(Xs>=0)){
    stop("Xs must contain positive integers")
  }
  if (!inherits(Xb, "matrix")){
    stop("Xb must be of class \"matrix\"")
  }
  if(!all(Xb==round(Xb)) | !is.numeric(Xb) | !all(Xb>=0)){
    stop("Xb must contain positive integers")
  }
  if (!(all(rownames(Rs) == rownames(Rb)))) {
    stop("Rownames for Rs and Rb are not identical")
  }
  if (!(all(rownames(Xs) == rownames(Xb)))) {
    stop("Rownames for Xs and Xb are not identical")
  }
  if (!(all(colnames(Rs) == colnames(Xs)))) {
    stop("Colnames for Rs and Xs are not identical")
  }
  if (!(all(colnames(Rb) == colnames(Xb)))) {
    stop("Colnames for Rb and Xb are not identical")
  }
  if (!all(alpha > 0) | !is.numeric(alpha)){
    stop("All values in vector alpha must be greater than 0")
  }
  if (length(alpha)!=nrow(Rs)){
    stop("alpha must be of length equal to the number of rows in Rs (number of mutations)")
  }
  if (!all(beta > 0) | !is.numeric(beta)){
    stop("All values in vector beta must be greater than 0")
  }
  if (length(beta)!=nrow(Rs)){
    stop("beta must be of length equal to the number of rows in Rs (number of mutations)")
  }
  if (length(kappa)!=1){
    stop("kappa must be of length 1")
  }
  if (kappa <= 0 | !is.numeric(kappa)){
    stop("kappa must be greater than 0")
  }
  if (length(tau)!=1){
    stop("tau must be of length 1")
  }
  if (tau <= 0 | !is.numeric(tau)){
    stop("tau must be greater than 0")
  }
  if (!all(Klist >= 3) | !is.numeric(Klist) | all(Klist != floor(Klist)) ){
    stop("All values in vector Klist must be positive integers with value greater than or equal to 3")
  }
  if (length(niter)!=1){
    stop("niter must be of length 1")
  }
  if (niter <= 0 | !is.numeric(niter) | niter!=round(niter)){
    stop("niter must be a positive integer")
  }
  if (length(nchains)!=1){
    stop("nchains must be of length 1")
  }
  if (nchains <= 0 | !is.numeric(nchains) | nchains!=round(nchains)){
    stop("nchains must be a positive integer")
  }
  if (length(thin)!=1){
    stop("thin must be of length 1")
  }
  if (thin <= 0 | !is.numeric(thin) | thin!=round(thin)){
    stop("thin must be a positive integer")
  }
  if (length(pburn)!=1){
    stop("pburn must be of length 1")
  }
  if (pburn<0 | pburn>=1 | !is.numeric(pburn)){
    stop("pburn must be between 0 (inclusive) and 1 (exclusive)")
  }
  if (!(length(ncores)%in%c(0,1))){
    stop("ncores must be of length 1")
  }
  if (!is.null(ncores)){
    if(ncores <= 0 | !is.numeric(ncores) | ncores!=round(ncores)){
      stop("ncores must be a positive integer")
    }
  }
  if(is.null(ncores)){
    ncores<-1/2*parallelly::availableCores()
  }
  if(ncores>parallelly::availableCores()){
    ncores<-parallelly::availableCores()
  }
  if(ncores==1){
    print("Executing code in serial using one core. If parallelization is desired, see argument 'ncores'")
  }else{
    print(paste("Executing code in parallel using", ncores, "cores."))
  }
  if (length(seed)!=1){
    stop("seed must be of length 1")
  }
  if (seed <= 0 | !is.numeric(seed) | seed!=round(seed)){
    stop("seed must be a positive integer")
  }

  # Number of iterations after thinning
  niter.thin<-niter/thin

  # Length of burn-in
  burn.len<-pburn*niter.thin

  # Number of iterations after burn-in
  niter.post.burn<-niter.thin-burn.len

  if(niter.post.burn<10){
    stop("There are less than 10 iterations remaining after thinning and removing the
         burn-in period. Please adjust your parameters accordingly.")
  }

  N <- ncol(Rs) # Number of cells
  S <- ncol(Rb) # Number of bulk samples
  M <- nrow(Rs) # Number of point mutations

  samples.out<-list()

  # Parallel processing across chains
  # Get operating system
  os<-canopy2:::get_os()

  # Get all possible combinations of number of chains and number of subclones
  grid.iter<-expand.grid(1:nchains,Klist)

  # If Windows
  # Note: Forking not supported in these cases, so must register cluster
  # https://www.oreilly.com/library/view/parallel-r/9781449317850/ch04.html
  if(os=="windows"){
    cl <- parallel::makePSOCKcluster(ncores)
    parallel::setDefaultCluster(cl)
    parallel::clusterExport(cl, list('initialsnv','getZ','logdBetaBinom','getPost',
                                     'nchains','grid.iter', 'niter', 'thin', 'niter.thin',
                                     'burn.len','Rs', 'Xs', 'Rb', 'Xb', 'S', 'N',
                                     'alpha', 'beta', 'kappa', 'tau', 'seed'),
                            envir = environment())
    parallel::clusterEvalQ(cl, {library(stats); library(ape); library(DirichletReg)})
    out.mcmc<-parallel::parLapply(cl, 1:nrow(grid.iter), function(x)
      canopy2:::MH_within_Gibbs(chain=grid.iter[x,1], K=grid.iter[x,2],
                                nchains=nchains, niter=niter, thin=thin,
                                niter.thin=niter.thin, burn.len=burn.len,
                                Rs=Rs, Xs=Xs, Rb=Rb, Xb=Xb, S=S, N=N,
                                alpha=alpha, beta=beta, kappa=kappa, tau=tau,
                                seed=seed))
    parallel::stopCluster(cl)
  }

  # If Darwin (macOS), Linux, or other Unix-based platform
  else{
    out.mcmc <- parallel::mclapply(1:nrow(grid.iter), function(x)
    {canopy2:::MH_within_Gibbs(chain=grid.iter[x,1], K=grid.iter[x,2],
                               nchains=nchains, niter=niter, thin=thin,
                               niter.thin=niter.thin, burn.len=burn.len,
                               Rs=Rs, Xs=Xs, Rb=Rb, Xb=Xb, S=S, N=N,
                               alpha=alpha, beta=beta, kappa=kappa, tau=tau,
                               seed=seed)}, mc.cores=ncores)
  }


  final.out<-list(samples=out.mcmc,
                  N=N,
                  S=S,
                  M=M,
                  Rs=Rs,
                  Rb=Rb,
                  Xs=Xs,
                  Xb=Xb,
                  alpha=alpha,
                  beta=beta,
                  kappa=kappa,
                  tau=tau,
                  Klist=Klist,
                  nchains=nchains,
                  pburn=pburn)
  class(final.out)<-"get_trees"
  return(final.out)

}
