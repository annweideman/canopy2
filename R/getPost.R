#' Compute the log posterior density
#'
#' Compute the log posterior density using the contributions from the
#' single-cell data (beta-binomial likelihoods) and bulk data (binomial
#' likelihood).
#'
#' @param tree an object of class \code{pyhlo} representing the inferred
#' phylogenetic tree.
#' @param Rb \eqn{M} (mutations) x \eqn{S} (bulk samples) matrix of bulk
#' alternative read counts. Rownames must match those of argument \code{Rs}.
#' @param Xb \eqn{M} (mutations) x \eqn{S} (bulk samples) matrix of bulk total
#' (benign + mutated) read counts. Rownames must match those of argument
#' \code{Xs}.
#' @param Rs \eqn{M} (mutations) x \eqn{N} (cells) matrix of single-cell
#' alternative read counts. Rownames must match those of argument \code{Rb}.
#' @param Xs \eqn{M} (mutations) x \eqn{N} (cells) matrix of single-cell
#' total (benign + mutated) read counts. Rownames must match those of argument
#' \code{Xb}.
#' @param alpha numeric vector of positive values with length \eqn{M}
#' (number of mutations) representing the gene activation rates.
#' @param beta numeric vector of positive values with length \eqn{M}
#' (number of mutations) representing the gene deactivation rates.
#' @param kappa a positive value used in the computation of the
#' sequencing error defined as \eqn{\kappa/(\kappa+\tau)}.
#' @param tau a positive value used in the computation of the sequencing error
#' defined as \eqn{\kappa/(\kappa+\tau)}.
#'
#' @return
#' A numeric representing the log posterior density.
#'
#' @examples
#' # Simulate read counts
#' sims.out<-simulate_data(N=15, S=2, M=5, alpha=0.1, beta=1.0,kappa=1, tau=999,
#'                         Ktrue=4, b.mindepth=30, b.maxdepth=50, sc.mindepth=80,
#'                         sc.maxdepth=120, scale=300, seed=8675309)
#'
#' # Estimate parameters for gene kinetics using the BPSC methodology
#' # Note: can also use get_burstiness_scale() for large datasets, but estimates
#' # are not as reliable
#' param.out<-get_burstiness_bpsc(counts=sims.out$G)
#'
#' # Subset the read counts to include only those mutations with estimable gene
#' # kinetics (all are estimable here, but included for utility when using real
#' # data)
#' Rs<-sims.out$Rs[param.out$id.g,] # single cell alternative read counts
#' Xs<-sims.out$Xs[param.out$id.g,] # single cell total read counts
#' Rb<-sims.out$Rb[param.out$id.g,] # bulk alternative read c
#'
#' N <- ncol(sims.out$Rs) # Number of cells
#' S <- ncol(sims.out$Rb) # Number of bulk samples
#' M <- nrow(sims.out$Rs) # Number of point mutations
#'
#' # Generate a tree
#' K<-4
#' text = paste0(paste0(paste(paste0("(", 1:(K - 1), ","),
#'                            collapse = ""), K), paste(rep(")", (K - 1)), collapse = ""), ";")
#' tree <- ape::read.tree(text = text); rm(text)
#'
#' # Generate point mutations (single-nucleotide varaints) along the tree branches
#' tree$snv<-initialsnv(tree, rownames(sims.out$Rs))
#'
#' # Get the Z matrix from tree and snv
#' tree$Z <- getZ(tree)
#'
#' # Generate the proportion matrices
#' # Bulk samples
#' Pb<-t(DirichletReg::rdirichlet(S, alpha=rep(1/K,K)))
#' rownames(Pb)<-paste0('clone',1:K)
#' colnames(Pb)<-colnames(sims.out$Rb)
#'
#' # Single cells
#' Ps<-stats::rmultinom(N, 1, prob = rep(1/K, K))
#' rownames(Ps)<-paste0('clone',1:K)
#' colnames(Ps)<-colnames(sims.out$Rs)
#' tree$Pb<-Pb; rm(Pb)
#' tree$Ps<-Ps; rm(Ps)
#'
#' # Evaluate the log posterior (with non-informative prior)
#' tree$Post<-getPost(tree, Rb=sims.out$Rb, Xb=sims.out$Xb, Rs=sims.out$Rs, Xs=sims.out$Xs,
#'                    alpha=param.out$alpha, beta=param.out$beta, kappa=1, tau=999)
#'
#' tree$Post
#'
#' @export

getPost<-function(tree, Rb, Xb, Rs, Xs, alpha, beta, kappa, tau){

  # Mixture of beta-binomial likelihoods for the single-cell data
  Qs=tree$Z%*%tree$Ps
  logPost=sum((logdBetaBinom(Rs, Xs, alpha, beta))*Qs+
              (logdBetaBinom(Rs, Xs, kappa, tau))*(1-Qs))
  
  # Combine with binomial likelihood (written up to a proportionality constant)
  # for the bulk data
  Qb=pmin(pmax(1/2*tree$Z%*%tree$Pb, 0.01),0.99)
  logPost=logPost+sum(lchoose(Xb,Rb) + Rb*log(Qb)+(Xb-Rb)*log(1-Qb))
  return(logPost)

}
