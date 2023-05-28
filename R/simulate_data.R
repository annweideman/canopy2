#' Simulate data for Canopy2
#'
#' Simulates data for Canopy2 from a beta-binomial distribution.
#'
#' @param N number of single cells
#' @param S number of bulk samples
#' @param M number of mutations (single nucleotide variants)
#' @param Ktrue true number of subclones
#' @param alpha a positive integer representing the gene activation rate. This
#' number will be replicated to form a vector of length \eqn{M} (number of mutations).
#' @param beta a positive integer representing the gene deactivation rate. This
#' number will be replicated to form a vector of length \eqn{M} (number of mutations).
#' @param kappa a positive value used in the computation of the sequencing error
#' defined as \eqn{\kappa/(\kappa+\tau)}. Defaults to 1.
#' @param tau a positive value used in the computation of the sequencing error
#' defined as \eqn{\kappa/(\kappa+\tau)}. Defaults to 999.
#' @param b.mindepth a positive value that represents the minimum sequencing
#' depth of the bulk data
#' @param b.maxdepth a positive value that represents the maximum sequencing
#' depth of the bulk data
#' @param sc.mindepth a positive value that represents the minimum sequencing
#' depth of the single-cell data
#' @param sc.maxdepth a positive value that represents the maximum sequencing
#' depth of the single-cell data
#' @param scale a positive value that represents the scaling factor for beta
#' distribution (the rate at which DNA transcribed into RNA)
#' @param seed a state (positive integer) to set the random number generation.
#' Defaults to 8675309.
#'
#' @return
#' A list containing: \code{true.tree}, the true/known phylogenetic tree,
#' \code{Rs}, a matrix of single-cell alternative read counts, \code{Rb}, a
#' matrix of bulk alternative read counts, \code{Xs}, a matrix of single-cell
#' total (benign + mutated) read counts, \code{Xb}, a matrix of bulk total
#' (benign + mutated) read counts, \code{alpha}, the mutation-specific gene
#' activation rates, \code{beta}, the mutation-specific gene deactivation rates,
#' and \code{G}, an \eqn{M} (mutation) x \eqn{N} (cell) matrix of single-cell gene
#' expression data.
#'
#' @export
#'
#' @examples
#' # Simulate read counts
#' sims.out<-simulate_data(N=30, S=4, M=15, alpha=0.1, beta=1, kappa=1, tau=999,
#'                         Ktrue=4, b.mindepth=30, b.maxdepth=50, sc.mindepth=80,
#'                         sc.maxdepth=120, scale=300, seed=8675309)
#'
#' # Run Canopy2 to get list of phylogenetic trees corresponding to all chains
#' # and all subclones
#' get.trees.out<-get_trees(Rs=sims.out$Rs, Rb=sims.out$Rb,
#'                          Xs=sims.out$Xs, Xb=sims.out$Xb,
#'                          alpha=sims.out$alpha, beta=sims.out$beta,
#'                          kappa=1, tau=999, Klist=3:5,
#'                          niter=10000, nchains=5, thin=10, pburn=0.1, seed=8675309)
#'
#' # Examine diagnostic plots
#' get_diagnostics(get.trees.out, project=NULL, outpath=NULL)
#'
#' # Get best tree across all chains and subclones via DIC
#' best.tree.out<-get_best_tree(get.trees.out)
#'
#' best.tree.out

simulate_data<-function(N, S, M, alpha, beta, kappa=1, tau=999, Ktrue,
                        b.mindepth, b.maxdepth, sc.mindepth, sc.maxdepth,
                        scale, seed=8675309){

  if (length(N)!=1){
    stop("N must be of length 1")
  }
  if (N <= 0 | !is.numeric(N)){
    stop("N must be a number greater than 0")
  }
  if (length(S)!=1){
    stop("S must be of length 1")
  }
  if (S <= 0 | !is.numeric(S)){
    stop("S must be a number greater than 0")
  }
  if (length(M)!=1){
    stop("M must be of length 1")
  }
  if (M <= 0 | !is.numeric(M)){
    stop("M must be a number greater than 0")
  }
  if (length(alpha)!=1){
    stop("alpha must be of length 1")
  }
  if (alpha <= 0 | !is.numeric(alpha)){
    stop("alpha must be a number greater than 0")
  }
  if (length(beta)!=1){
    stop("beta must be of length 1")
  }
  if (beta <= 0 | !is.numeric(beta)){
    stop("beta must be a number greater than 0")
  }
  if (length(kappa)!=1){
    stop("kappa must be of length 1")
  }
  if (kappa <= 0 | !is.numeric(kappa)){
    stop("kappa must be a number greater than 0")
  }
  if (length(tau)!=1){
    stop("tau must be of length 1")
  }
  if (tau <= 0 | !is.numeric(tau)){
    stop("tau must be a number greater than 0")
  }
  if (length(Ktrue)!=1){
    stop("Ktrue must be of length 1")
  }
  if (Ktrue <= 0 | !is.numeric(Ktrue)){
    stop("Ktrue must be a number greater than 0")
  }
  if (length(b.mindepth)!=1){
    stop("b.mindepth must be of length 1")
  }
  if (b.mindepth <= 0 | !is.numeric(b.mindepth)){
    stop("b.mindepth must be a number greater than 0")
  }
  if (length(b.maxdepth)!=1){
    stop("b.maxdepth must be of length 1")
  }
  if (b.maxdepth <= 0 | !is.numeric(b.maxdepth)){
    stop("b.maxdepth must be a number greater than 0")
  }
  if (b.mindepth>=b.maxdepth){
    stop("b.mindepth must be strictly less than b.maxdepth")
  }
  if (length(sc.mindepth)!=1){
    stop("sc.mindepth must be of length 1")
  }
  if (sc.mindepth <= 0 | !is.numeric(sc.mindepth)){
    stop("sc.mindepth must be a number greater than 0")
  }
  if (length(sc.maxdepth)!=1){
    stop("sc.maxdepth must be of length 1")
  }
  if (sc.maxdepth <= 0 | !is.numeric(sc.maxdepth)){
    stop("sc.maxdepth must be a number greater than 0")
  }
  if (sc.mindepth>=sc.maxdepth){
    stop("sc.mindepth must be strictly less than sc.maxdepth")
  }
  if (length(scale)!=1){
    stop("scale must be of length 1")
  }
  if (scale <= 0 | !is.numeric(scale)){
    stop("scale must be a number greater than 0")
  }
  if (length(seed)!=1){
    stop("seed must be of length 1")
  }
  if (seed <= 0 | !is.numeric(seed) | seed!=round(seed)){
    stop("seed must be a positive integer")
  }
  if (Ktrue>M){
    stop("Ktrue must be less than or equal to M")
  }
  if (Ktrue>N){
    stop("Ktrue must be less than or equal to N")
  }

  set.seed(seed)

  # Generate a tree
  text <- paste0(paste0(paste(paste0("(", 1:(Ktrue - 1), ","),
                             collapse = ""), Ktrue), paste(rep(")", (Ktrue - 1)),
                                                            collapse = ""), ";")
  tree <- ape::read.tree(text = text); rm(text)
  tree$edge # The edge of the tree

  tree$Z <- matrix(c(0,0),ncol=2) # initialize while loop

  # Check for duplicated columns in true Z and rerun sampling if present
  while(any(duplicated(t(tree$Z)))){

    # Generate point mutations (SNV: single-nucleotide variant)
    # along the tree branches
    tree$snv <- initialsnv(tree, paste0('snv',1:M))

    # Get the Z matrix from the tree and SNVs
    tree$Z <- getZ(tree)# M mutations x K clones

  }

  # Generate the proportion matrices
  # Bulk samples
  Pb<-t(DirichletReg::rdirichlet(S, alpha=rep(1/Ktrue,Ktrue)))
  rownames(Pb)<-paste0('clone',1:Ktrue)
  colnames(Pb)<-paste0('bulk',1:S)
  apply(Pb,2,sum)

  # Single cells
  Ps <- stats::rmultinom(N, 1, prob = rep(1/Ktrue, Ktrue))
  rownames(Ps) <- paste0('clone',1:Ktrue)
  colnames(Ps) <- paste0('cell',1:N)
  apply(Ps,2,sum)
  tree$Pb <- Pb; rm(Pb)
  tree$Ps <- Ps; rm(Ps)

  # Generate the input matrices Rb, Xb from a simple binomial model
  Qb <- tree$Z%*%tree$Pb
  Rb <- matrix(nrow=M, ncol=S)
  rownames(Rb) <- paste0('snv',1:M)
  colnames(Rb) <- paste0('bulk',1:S)
  Xb <- Rb
  # Total coverage between mindepth - maxdepth
  Xb[1:length(Xb)] <- round(stats::runif(n = length(Xb), min=b.mindepth, max=b.maxdepth))
  true.prob <- pmax(1/2*Qb,0.01)
  for(i in 1:nrow(Rb)){
    for(t in 1:ncol(Rb)){
      Rb[i,t] <- stats::rbinom(1, size = Xb[i,t], prob = true.prob[i,t])
    }
  }

  # Generate the input matrices Rs, Xs from a beta-binomial model
  Qs <- tree$Z%*%tree$Ps
  true.prob <- pmax(Qs, 0.01)

  eps <- stats::rbeta(M, kappa, tau)
  alpha <- rep(alpha, M)
  beta <- rep(beta, M)
  scale <- rep(scale, M)
  names(eps) <- names(alpha)<-names(beta)<-names(scale)<-paste0('snv',1:M)

  Rs <- matrix(nrow=M, ncol=N)
  rownames(Rs) <- paste0('snv',1:M)
  colnames(Rs) <- paste0('cell',1:N)
  Xs <- Rs
  # Total coverage between mindepth - maxdepth
  Xs[1:length(Xs)] <- round(stats::runif(n = length(Xs), min=sc.mindepth, max=sc.maxdepth))

  for(i in 1:nrow(Rs)){
    for(j in 1:ncol(Rs)){
      if(Qs[i,j]==1){
        Gamma.ij <- stats::rbeta(1, alpha[i], beta[i])
        Rs[i,j] <- stats::rbinom(1, size = Xs[i,j], prob = Gamma.ij)
      } else{
        Rs[i,j] <- stats::rbinom(1, size = Xs[i,j], prob = eps[i])
      }
    }
  }

  # Sort Ps to aggregate 1's for better visualization
  temp.mat<-matrix(NA,nrow=nrow(tree$Ps),ncol=1)
  rownames(temp.mat)<-rownames(tree$Ps)
  id.ones<-lapply(1:nrow(tree$Ps),function(x) which(tree$Ps[x,]==1))

  # Generate mini binary matrices that are concatenated to create final matrix
  for(list.id in length(id.ones):1){
    temp.mat0<-matrix(0,nrow=nrow(tree$Ps),ncol=length(id.ones[[list.id]]))
    temp.mat0[list.id,]<-1
    colnames(temp.mat0)<-names(id.ones[[list.id]])
    temp.mat<-cbind(temp.mat,temp.mat0)
  }

  # Remove vector of NAs
  temp.mat<-temp.mat[,-1]
  # Replace Ps in final tree with new Ps that has 1's along main diagonal
  tree$Ps<-temp.mat

  # Generate random variates from a beta-poisson to simulate gene
  # expression data
  G <- t(sapply(1:length(alpha), function(i)
           rBetaPois(n=N,alpha=alpha[i],beta=beta[i],scale=scale[i])))
  rownames(G) <- paste0("snv", 1:nrow(G))
  colnames(G) <- paste0("cell", 1:ncol(G))

  return(list(true.tree=tree, Rs=Rs, Rb=Rb, Xs=Xs, Xb=Xb, alpha=alpha, beta=beta, G=G))

}