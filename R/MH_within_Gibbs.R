#' Internal function
#'
#' Gibbs with nested Metropolis-Hastings. See Alg. 2 in the main text for
#' further details.
#'
#' @examples \dontrun{
#' MH_within_Gibbs(chain=1, K=4, seed=8675309)
#' }
#' @keywords internal
#'
MH_within_Gibbs<-function(chain, K, nchains, niter, thin, niter.thin, burn.len,
                          Rs, Xs, Rb, Xb, S, N, alpha, tau, seed=8675309){

  # Print message to console
  print(paste("Running chain", chain, "out of", nchains, "..."))

  tree <- canopy2:::initialize_tree(seedling = (chain + seed),
                                    K=K,
                                    Rs=Rs,
                                    S=S,
                                    Rb=Rb,
                                    N=N,
                                    Xb=Xb,
                                    Xs=Xs,
                                    alpha=alpha,
                                    tau=tau)
  tree.list <- list()
  accept <- numeric() # Initialize an empty numeric vector to store acceptance values

  for (iter in 1:niter) {

    tosample <- sample.int(1, n = 3)

    if (tosample == 1) {
      sampZout <- canopy2:::sampZ(tree,Rb, Xb, Rs, Xs, alpha, beta, kappa, tau)
      tree <- sampZout[[1]]
      accept <- c(accept, sampZout[[2]])
    } else if (tosample == 2) {
      sampPbout <- canopy2:::sampPb(tree,Rb, Xb, Rs, Xs, alpha, beta, kappa, tau)
      tree <- sampPbout[[1]]
      accept <- c(accept, sampPbout[[2]])
    } else {
      sampPsout <- canopy2:::sampPs(tree,Rb, Xb, Rs, Xs, alpha, beta, kappa, tau)
      tree <- sampPsout[[1]]
      accept <- c(accept, sampPsout[[2]])
    }

    # Store only thinned samples from posterior
    if (iter %% thin == 0) {
      tree.list[[length(tree.list) + 1]] <- tree
    }
  }

  # Grab posterior probabilities
  post <- sapply(1:niter.thin, function(x) tree.list[[x]]$Post)

  # Remove burn-in from samples
  tree.list <- tree.list[-c(1:burn.len)]
  post <- post[-c(1:burn.len)]
  accept <- accept[-c(1:burn.len)]
  accept.rate <- sum(accept) / length(accept)

  return(list("K" = K, "tree" = tree.list, "posteriors" = post, "acceptance rate" = accept.rate))

}
