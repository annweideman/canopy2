#' Sample space of phylogenetic trees
#'
#' Produce a list of phylogenetic trees in a sample space of size
#' number of chains x number of subclones. Sampling is performed using
#' Metropolis-within-Gibbs. Output is coded as class \code{get_trees}.
#'
#' @param Rs mutation x cell matrix of integers. Denotes the number 
#' of alternative allele reads corresponding to variant m in cell n. 
#' Rownames must match those of argument \code{Rb}
#' @param Rb mutation x bulk sample matrix of integers. Denotes the number 
#' of alternative allele reads corresponding to variant m in bulk sample s. 
#' Rownames must match those of argument \code{Rs}
#' @param Xs mutation x cell matrix of integers. Denotes the number 
#' of total reads corresponding to variant m in cell n. 
#' Rownames must match those of argument \code{Xb}
#' @param Xb mutation x bulk sample matrix of integers. Denotes the number 
#' of total reads corresponding to variant m in bulk sample s. 
#' Rownames must match those of argument \code{Xs}
#' @param alpha numeric vector of values greater than 0 with length \eqn{M}
#' (number of mutations) representing the gene activation rate. These can be
#' estimated using the 
#' @param beta numeric vector of values greater than 0 with length \eqn{M}
#' (number of mutations) representing the gene deactivation rate
#' @param kappa a positive value used in the computation of the
#' sequencing error defined as \eqn{\kappa/(\kappa+\tau)}
#' @param tau a positive value used in the computation of the sequencing error
#' defined as \eqn{\kappa/(\kappa+\tau)}
#' @param Klist numeric vector containing the possible numbers of subclones
#' @param niter number of iterations of MCMC. Defaults to 10000.
#' @param nchains number of chains for MCMC. Defaults to 20.
#' @param thin a number representing the increment at which to store MCMC output.
#'  For example, if \code{thin=10} then every 10th iteration is stored. Defaults
#'  to 10.
#' @param pburn a decimal denoting the percentage of burn-in to store. Defaults
#' to 0.10 (10%). 
#' @param seed a state (positive integer) to set the random number generation.
#' Defaults to 1.
#'
#' @return
#' @export
#'
#' @examples
get_trees<-function(Rs, Rb, Xs, Xb, alpha, beta, kappa, tau,
                    Klist, niter=10000, nchains=20, thin=10, pburn=0.1, seed=1){

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
    stop("alpha must be of length equal to the number of rows in Rs (number of
         mutations)")
  }
  if (!all(beta > 0) | !is.numeric(beta)){
    stop("All values in vector beta must be greater than 0")
  }
  if (length(beta)!=nrow(Rs)){
    stop("beta must be of length equal to the number of rows in Rs (number of
         mutations)")
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

  # Loop over number of subclones
  for(K in Klist){

    # Print message to console
    print(paste("Sample in tree space with", K, "subclones"))

   #----------------------------------------------------------------------------
   # Estimating unknowns
   #----------------------------------------------------------------------------

    # The unknowns are Pb, Ps, Z

    #-----------------------------------------------
    # Initialization
    #-----------------------------------------------
    initialize<-function(seedling){
      
      set.seed(seedling)

      # Generate a tree
      text = paste0(paste0(paste(paste0("(", 1:(K - 1), ","),
                                 collapse = ""), K), paste(rep(")", (K - 1)), collapse = ""), ";")
      tree <- ape::read.tree(text = text); rm(text)
      tree$edge # The edge of the tree

      # Generate point mutations (SNV: single-nucleotide alteration) along the tree branches
      tree$snv=initialsnv(tree, rownames(Rs))

      # Get the Z matrix from tree and snv
      tree$Z = getZ(tree)
      tree$Z # M mutations x K clones

      # Generate the proportion matrix
      # Bulk samples
      Pb=t(DirichletReg::rdirichlet(S, alpha=rep(1/K,K)))
      rownames(Pb)=paste0('clone',1:K)
      colnames(Pb)=colnames(Rb)
      apply(Pb,2,sum)
      
      # Single cells
      Ps=stats::rmultinom(N, 1, prob = rep(1/K, K))
      rownames(Ps)=paste0('clone',1:K)
      colnames(Ps)=colnames(Rs)
      apply(Ps,2,sum)
      tree$Pb=Pb; rm(Pb)
      tree$Ps=Ps; rm(Ps)

      # Evaluate the posterior (with non-informative prior)
      tree$Post=getPost(tree, Rb, Xb, Rs, Xs, alpha, beta, kappa, tau)

      return(tree)
    }

    #---------------------------------------------------
    # Metropolis sampling (see Algorithm 1, main text)
    #---------------------------------------------------
    sampZ=function(tree,...){
      M = nrow(tree$snv) # number of mutations
       for(i in 1:M){ # Update each mutation through a loop
        snv.new = tree$snv
        # sample from edges excluding the leftmost and the current edge
        # the probability for each edge to be selected is the same (i.e., symmetric)
        snv.edge.i = sample(setdiff(2:nrow(tree$edge), tree$snv[i,4]), size = 1)
        snv.new[i, 2:4] = c(tree$edge[snv.edge.i, 1:2], snv.edge.i)
        tree.new=tree
        tree.new$snv=snv.new
        tree.new$Z=getZ(tree.new)

        tree.new$Post=getPost(tree.new, Rb, Xb, Rs, Xs, alpha, beta, kappa, tau)
        r=exp(tree.new$Post-tree$Post)
        if(r>= stats::runif(1)){
          tree=tree.new
          accept<-c(accept,1)
        } else{
          tree=tree
          accept<-c(accept,0)
        }
      }
      return(list(tree,accept))
    }

    sampPs=function(tree,...){
      N=ncol(tree$Ps) # number of cells
      K=nrow(tree$Ps) # number of clones
      for(j in 1:N){ # Update each cell through a loop
          Ps.new=tree$Ps
          # sample from clones excluding the current clone
          # the probability for each clone to be selected is the same (i.e., symmetric)
          Ps.new[,j]=stats::rmultinom(1,1,prob=(1-tree$Ps[,j])/(K-1))
          tree.new=tree
          tree.new$Ps=Ps.new
          tree.new$Post=getPost(tree.new, Rb, Xb, Rs, Xs, alpha, beta, kappa, tau)
          
          r=exp(tree.new$Post-tree$Post)
          if(r>= stats::runif(1)){
            tree=tree.new
            accept<-c(accept,1)
          } else{
            tree=tree
            accept<-c(accept,0)
          }
        }
      return(list(tree,accept))
    }
    
    sampPb=function(tree,...){
      S=ncol(tree$Pb) # number of bulk samples
      K=nrow(tree$Pb) # number of clones
      for(t in 1:S){ # Update each bulk sample through a loop
          Pb.new=tree$Pb
          # dirichlet distribution with concentration parameter equal to the current proportion
          Pb.new[,t]= gtools::rdirichlet(1,alpha=Pb.new[,t])
          tree.new=tree
          tree.new$Pb=Pb.new
          tree.new$Post=getPost(tree.new, Rb, Xb, Rs, Xs, alpha, beta, kappa, tau)
          r=exp(tree.new$Post+log(gtools::ddirichlet(tree$Pb[,t],alpha=Pb.new[,t]))
                -tree$Post-log(gtools::ddirichlet(Pb.new[,t],alpha=tree$Pb[,t])))
          if(r>= stats::runif(1)){
            tree=tree.new
            accept<-c(accept,1)
          }else{
            tree=tree
            accept<-c(accept,0)
          }
        #}
      }
      return(list(tree,accept))
    }


    accept<-c()

    #--------------------------------------------------------------------
    # Gibbs with nested Metropolis-Hastings (see Algorithm 2, main text)
    #--------------------------------------------------------------------

    out.mcmc<-lapply(1:nchains, function(chain){

      # Print message to console
      print(paste("Running chain",chain, "out of", nchains,"..."))

      tree<-initialize(seedling=chain)
      tree.list<-list()

      for(iter in 1:niter){

        tosample=sample.int(1, n=3)
        
        if(tosample==1){
          sampZout=sampZ(tree)
          tree=sampZout[[1]]
          accept=c(accept,sampZout[[2]])
        }else if(tosample==2){
          sampPbout=sampPb(tree)
          tree=sampPbout[[1]]
          accept=c(accept,sampPbout[[2]])
        }else{
          sampPsout=sampPs(tree)
          tree=sampPsout[[1]]
          accept=c(accept,sampPsout[[2]])
        }

        # Store only thinned samples from posterior
        if(iter %% thin==0){
          tree.list[[length(tree.list)+1]]<-tree
        }
      }

      # Grab posterior probabilities
      post<-sapply(1:niter.thin, function(x) tree.list[[x]]$Post)

      # Remove burn-in from samples
      tree.list<-tree.list[-c(1:burn.len)]
      post<-post[-c(1:burn.len)]
      accept<-accept[-c(1:burn.len)]
      accept.rate<-sum(accept)/length(accept)
      
      # Sort Ps to aggregate 1's for better visualization
      temp.mat<-matrix(NA,nrow=nrow(tree$Ps),ncol=1)
      rownames(temp.mat)<-rownames(tree$Ps)
      id.ones<-lapply(1:nrow(tree$Ps),function(x) which(tree$Ps[x,]==1))
      
      # Generate mini binary matrices that are concatenated to create final matrix
      for(list.id in 1:length(id.ones)){
        temp.mat0<-matrix(0,nrow=nrow(tree$Ps),ncol=length(id.ones[[list.id]]))
        temp.mat0[list.id,]<-1
        colnames(temp.mat0)<-names(id.ones[[list.id]])
        temp.mat<-cbind(temp.mat,temp.mat0)
      }
      
      # Remove vector of NAs
      temp.mat<-temp.mat[,-1]
      
      return(list("K"=K,"tree"=tree.list,"posteriors"=post,"acceptance rate"=accept.rate))
      }
      )

     samples.out<-c(samples.out,out.mcmc)
    
  }

  final.out<-list(samples=samples.out,
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
