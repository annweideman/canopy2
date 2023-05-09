#-------------------------------------------------------------------------------
# Script to store helper functions
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Function for computing distance metric for absolute reconstruction error
#-------------------------------------------------------------------------------

dist_metric_fun <- function(truth, inferred){

  # Find the maximum number of columns
  max.ncol <- max(ncol(truth),ncol(inferred))

  # Store all possible permutations of column numbers
  perms <- gtools::permutations(n=max.ncol,r=max.ncol,v=1:max.ncol)
  # If the matrices differ in size, then remove the overage from the
  # permutations
  perms <- perms[,1:min(ncol(inferred),ncol(truth))]

  # Generate all columnwise permutations of the matrix with the most number
  # of columns
  # If number of columns in the inference exceeds or equals the number of
  # columns in the truth
  if(ncol(truth) <= ncol(inferred)){
    mat.perms <- lapply(1:nrow(perms), function(x) inferred[,perms[x,]])
    # Compute the distance between the truth and each permuted inference
    mat.dist <- lapply(1:length(mat.perms), function(x)
      sum(colSums(abs(truth-mat.perms[[x]]))))
    # If number of columns in the truth exceeds number of columns in the inference
  }else{
    mat.perms <- lapply(1:nrow(perms), function(x) truth[,perms[x,]])
    # Compute the distance between each permuted truth and the inference
    mat.dist <- lapply(1:length(mat.perms), function(x)
      sum(colSums(abs(mat.perms[[x]]-inferred))))
  }

  # Compute the minimum distance metric
  dist.metric <- min(unlist(mat.dist))

  # Store the permutation corresponding to the minimum distance metric
  mat.perm.min <- mat.perms[[which.min(unlist(mat.dist))]]

  return(list(dist.metric=dist.metric,mat.perm.min=mat.perm.min))
}

#-------------------------------------------------------------------------------
# Function to compute absolute reconstruction error
#-------------------------------------------------------------------------------
are_fun <- function(true.tree, inferred.tree){

  N<-ncol(inferred.tree$Ps)
  #--------------------------------------------------------------
  # Transpose when needed to have subclones as columns
  #--------------------------------------------------------------
  true.Z <- true.tree$Z; inferred.Z <- inferred.tree$Z
  true.Ps <- t(true.tree$Ps); inferred.Ps <- t(inferred.tree$Ps)
  true.Pb <- t(true.tree$Pb); inferred.Pb <- t(inferred.tree$Pb)

  #-------------------------------------------
  # Compute distance metric for each component
  #-------------------------------------------
  dist.Z <- dist_metric_fun(truth=true.Z, inferred=inferred.Z)$dist.metric
  dist.Ps <- dist_metric_fun(truth=true.Ps, inferred=inferred.Ps)$dist.metric
  out.dist.Pb <- dist_metric_fun(truth=true.Pb, inferred=inferred.Pb)
  dist.Pb <- out.dist.Pb$dist.metric
  perm.min.Pb <- out.dist.Pb$mat.perm.min

  #--------------------------------------------------------------
  # Compute extra error from additional column(s), if applicable
  #--------------------------------------------------------------
  ### Extra error for Z
  ncol.diff.Z <- abs(ncol(inferred.Z)-ncol(true.Z))

  # If number of columns in the inference differs from number of columns in
  # the truth
  if(ncol.diff.Z!=0){
    # Add additional error that is equal to the length of the extra
    # column(s), M, times the number of extra column(s)
    e.Z <- M*ncol.diff.Z
    # Otherwise, set error due to extra columns to 0
  }else{e.Z <- 0}

  ### Extra error for Ps
  ncol.diff.Ps <- abs(ncol(inferred.Ps)-ncol(true.Ps))

  # If number of columns in the inference differs from number of columns in
  # the truth
  if(ncol.diff.Ps!=0){
    # Add additional error that is equal to the length of the extra
    # column(s), N, times the number of extra column(s)
    e.Ps <- N*ncol.diff.Ps
    # Otherwise, set error due to extra columns to 0
  }else{e.Ps <- 0}

  ### Extra error for Pb
  ncol.diff.Pb <- abs(ncol(inferred.Pb)-ncol(true.Pb))

  # If number of columns in the inference differs from number of columns in
  # the truth
  if(ncol.diff.Pb!=0){

    # Add additional error that is equal to the sum of the elements in the extra
    # column(s)
    if(ncol(inferred.Pb)>ncol(true.Pb)){
      e.Pb <- sum(inferred.Pb[,!(data.frame(inferred.Pb) %in% data.frame(perm.min.Pb))])
    }else{
      e.Pb <- sum(perm.min.Pb[,!(data.frame(perm.min.Pb) %in% data.frame(inferred.Pb))])
    }

    # Otherwise, set error due to extra columns to 0
  }else{e.Pb <- 0}

  #-------------------------------------------
  # Estimate the ARE for each component
  #-------------------------------------------

  are.Z <- dist.Z+e.Z
  are.Ps <- dist.Ps+e.Ps
  are.Pb <- dist.Pb+e.Pb

  #-------------------------------------------
  # Convert the estimated AREs to percentages
  #-------------------------------------------

  are.Z <- are.Z/max(prod(dim(inferred.Z)),prod(dim(true.Z)))*100
  are.Ps <- are.Ps/max(prod(dim(inferred.Ps)),prod(dim(true.Ps)))*100
  are.Pb <- are.Pb/(2*max(nrow(inferred.Pb), nrow(true.Pb)))*100

  return(list(are.Z=are.Z, are.Ps=are.Ps, are.Pb=are.Pb))

}

#-------------------------------------------------------------------------------
# Function to compute absolute reconstruction error for Z
#-------------------------------------------------------------------------------

are_z_fun <- function(true.Z, inferred.Z){

  # Check for duplicated columns in inferred Z and remove if present ONLY
  # if inference is larger than truth (i.e., number of columns in inference
  # exceeds number of columns in truth)
  if(any(duplicated(t(inferred.Z))) & ncol(inferred.Z) > ncol(true.Z)){
    # Which columns are duplicated
    id.Z.dup <- which(duplicated(t(inferred.Z)))
    # Remove duplicates
    inferred.Z <- inferred.Z[,-id.Z.dup]
  }

  M <- nrow(inferred.Z)
  K <- ncol(inferred.Z)

  #--------------------------------
  # Compute distance metric for Z
  #--------------------------------
  dist.Z <- dist_metric_fun(truth=true.Z, inferred=inferred.Z)$dist.metric

  #--------------------------------------------------------------
  # Compute extra error from additional column(s), if applicable
  #--------------------------------------------------------------
  # Extra error for Z
  ncol.diff.Z <- abs(ncol(inferred.Z)-ncol(true.Z))

  # If number of columns in the inference differs from number of columns in
  # the truth
  if(ncol.diff.Z!=0){
    # Add additional error that is equal to the length of the extra
    # column(s), M, times the number of extra column(s)
    e.Z <- M*ncol.diff.Z
    # Otherwise, set error due to extra columns to 0
  }else{e.Z <- 0}

  #-------------------------
  # Estimate the ARE for Z
  #-------------------------

  are.Z <- dist.Z+e.Z

  #--------------------------------------------------
  # Convert the estimated ARE for Z to a percentage
  #--------------------------------------------------

  # Divide by the total number of positions in the largest matrix
  are.Z <- are.Z/max(prod(dim(true.Z)),prod(dim(inferred.Z)))*100

  are.Z
}

#-------------------------------------------------------------------------------
# Compute deviance information criterion (DIC) according to Spiegelhalter defn
#-------------------------------------------------------------------------------
dic_fun <- function(logPost.mean, Z.mean, Ps.mean, Pb.mean, Rb,
                    Xb, Rs, Xs, alpha, beta, kappa, tau) {

  # Mixture of beta-binomial likelihoods for the single cell data
  Qs <- Z.mean%*%Ps.mean
  logPost <- sum((logdBetaBinom(Rs, Xs, alpha, beta))*Qs+
                   (logdBetaBinom(Rs, Xs, kappa, tau))*(1-Qs))

  # Combine with binomial likelihood (written up to a proportionality constant)
  # for the bulk data
  Qb <- pmin(pmax(1/2*Z.mean%*%Pb.mean, 0.01),0.99)
  logPost <- logPost+sum(Rb*log(Qb)+(Xb-Rb)*log(1-Qb))

  #logPost.mean <- mean(sapply(1:length(samples$tree), function(x)
  #                     getPost(samples$tree[[x]], Rb, Xb, Rs, Xs,
  #                             alpha, beta, kappa, tau)))

  # Effective number of parameters
  pD <- -2*logPost.mean-(-2*logPost)
  # Spiegelhalter defn of DIC
  DIC <- -2*logPost+2*pD

  return(DIC)
}

#-------------------------------------------------------------------------------
# Flatten a list of lists
#-------------------------------------------------------------------------------
flatten_list <- function(x){
  morelists <- sapply(x, function(xprime) class(xprime)[1]=="list")
  out <- c(x[!morelists], unlist(x[morelists], recursive=FALSE))
  if(sum(morelists)){
    Recall(out)
  }else{
    return(out)
  }
}

#-------------------------------------------------------------------------------
# Function to further process gene expression data:
# 1. Map Ensembl IDs to gene symbols
# 2. Remove any duplicates for Ensembl ID
#-------------------------------------------------------------------------------

get_gene_expression <- function(featurecounts, build){

  # Map Ensembl IDs to gene symbols
  ensembl <- biomaRt::useEnsembl(biomart="ensembl",
                                 dataset="hsapiens_gene_ensembl",
                                 GRCh=build)
  gene.names <- rownames(featurecounts)
  bm.out <- biomaRt::getBM(filters= "ensembl_gene_id",
                           attributes= c("ensembl_gene_id","hgnc_symbol"),
                           values=gene.names, mart=ensembl)

  # Remove any duplicates for Ensembl ID
  # Note: Ensembl arbitrarily picks a HGNC synonym for the summary if repeat
  # Ensembl IDs
  bm.out <- bm.out[!duplicated(bm.out$ensembl_gene_id), ]
  rownames(bm.out) <- bm.out$ensembl_gene_id

  # Merge with gene expression data
  featurecounts <- merge(bm.out,featurecounts,by="row.names",all=T)[,-1]

  featurecounts

}

#-------------------------------------------------------------------------------
# Compute the posterior density
#-------------------------------------------------------------------------------
getPost=function(tree, Rb, Xb, Rs, Xs, alpha, beta, kappa, tau){

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

#-------------------------------------------------------------------------------
# Generate binary clonal tree configuration matrix, Z
#-------------------------------------------------------------------------------
getZ <- function(tree) {

  K <- (nrow(tree$edge) + 2)/2 #Number of clones
  M <- nrow(tree$snv) #Number of mutations

  # Create a matrix of zeros
  Z <- matrix(nrow = M, ncol = K, data = 0)
  rownames(Z) <- rownames(tree$snv)
  colnames(Z) <- paste("clone", 1:K, sep = "")

  clonal.snv <- vector("list", K)
  for (tip in 2:K) {

    # Set child node to be the tip (clone) of interest
    child.node <- tip
    # Find parent node corresponding to child node
    parent.node <- tree$edge[which(tree$edge[, 2] == child.node), 1]

    # As long as the parent node has value greater than the number of clones
    while (parent.node > K) {

      # Find SNVs where parent and child node correspond to start and end node
      snvtemp <- intersect(which(tree$snv[, 2] == parent.node),
                           which(tree$snv[, 3] == child.node))

      # If an SNV exists along that branch, store in appropriate list
      if (length(snvtemp) > 0) {
        clonal.snv[[tip]] <- c(clonal.snv[[tip]], snvtemp)
      }

      # Set child node equal to parent node
      child.node <- parent.node
      # If child node has reached the top of the tree, then end while loop
      if (child.node == (K + 1))
        break

      # Find parent node corresponding to updated child node
      parent.node <- tree$edge[which(tree$edge[, 2] == child.node), 1]
    }
  }
  clonal.snv[[1]] <- 0
  for (ki in 2:K) {
    # Impute value of 1 for positions where SNVs exist
    Z[clonal.snv[[ki]],ki] <- 1
  }
  return(Z)
}

#-------------------------------------------------------------------------------
# Generate point mutations (SNVs: single-nucleotide variants) along the tree
# branches
#-------------------------------------------------------------------------------
initialsnv <- function (tree, snv.name) {
  # Number of point mutations
  snv.no <- length(snv.name)
  # If tree has more than two edges, then randomly place SNVs along these edges
  if (nrow(tree$edge) > 2) {
    snv.edge <- sample(2:nrow(tree$edge), size = snv.no, replace = T)
    # Otherwise, put SNVs on the only other edge besides the edge corresponding to
    # the normal clone
  } else {
    snv.edge <- rep(2, snv.no)
  }
  snv.mat <- cbind(1:snv.no, tree$edge[snv.edge, ], snv.edge)
  colnames(snv.mat) <- c("snv", "snv.strt.node", "snv.end.node", "snv.edge")
  rownames(snv.mat) <- snv.name
  return(snv.mat)
}

#-------------------------------------------------------------------------------
# Compute the log beta-binomial density up to a proportionality constant
#-------------------------------------------------------------------------------
logdBetaBinom=function(r, x, shape1, shape2){
  lchoose(x,r)+lbeta(r+shape1, x-r+shape2)-lbeta(shape1, shape2)
}

#-------------------------------------------------------------------------------
# Plot the clonal tree configuration matrix
#-------------------------------------------------------------------------------
plot_Z <- function(tree) {
  node_total <- max(tree$edge)
  node_shown <- ncol(tree$Z)
  node_hidden <- node_total - node_shown
  rowsums_Ps <- rowSums(tree$Ps)
  tree$tip.label[seq_len(node_shown)] <- paste0(
    "C", seq_len(node_shown),
    ": ", round(rowsums_Ps/sum(rowsums_Ps) * 100, digits = 0), "%")
  mut_id_all <- tree$Z %*% (2**seq(ncol(tree$Z), 1))
  mut_id_all <- seq(length(unique(mut_id_all)), 1)[as.factor(mut_id_all)]

  count_mut_ids<-0
  branch_ids <- rep(NA,node_total)

  # alternate binding of vectors
  trunk_ids<-(node_shown+1):node_total
  stem_ids<-1:node_shown
  alt_inds<-c(trunk_ids, stem_ids)[order(c(seq_along(trunk_ids)*2 - 1, seq_along(stem_ids)*2))]

  for (i in alt_inds) {
    mut_num <- sum(tree$snv[, 3] == i)
    if (mut_num == 0) {
      if (i == node_shown + 1) {
        branch_ids[i] <- "Root"
      } else {
        branch_ids[i] <- ""
      }
    } else{
      count_mut_ids<-count_mut_ids+1
      branch_ids[i] <-
        paste0("M", count_mut_ids, ": ", mut_num, " SNVs")
    }
  }

  pt <- ggtree::ggtree(tree)
  pt <- pt + ggplot2::geom_label(ggplot2::aes_string(x = "branch"),
                                 label = branch_ids, color = "firebrick"
  )

  pt <- pt + ggtree::geom_tiplab(hjust = 0.39, vjust = 0.85) +
    ggplot2::scale_x_reverse() + ggplot2::coord_flip()

  #top, right, bottom, left
  pt+ggplot2::theme(plot.margin = grid::unit(c(0,2,0,0), "cm"))
}

#-------------------------------------------------------------------------------
# Plot the phylogenetic tree with Z, Ps, and Pb
#-------------------------------------------------------------------------------

plot_tree <- function(tree, annovar=NULL, save.muts=F, save.plot=F,
                      outpath=NULL, filename=NULL,...){

  snvedge <- rep(NA, nrow(tree$snv))
  for (k in 1:nrow(tree$snv)) {
    snvedge[k] <- intersect(which(tree$edge[, 1] == tree$snv[k, 2]),
                            which(tree$edge[, 2] == tree$snv[k, 3]))
  }

  edge.label <- sort(unique(snvedge))
  mut.list=matrix(nrow=length(edge.label),ncol=1)
  if(is.null(annovar)){snv.name=rownames(tree$Z)
  }else{
    snv.name=rownames(Rs)}

  # Generate list of mutation clusters
  for (i in 1:length(edge.label)) {
    gene <- snv.name[which(snvedge == edge.label[i])]
    mut.list[i,1]=paste0("M", i, ": ", paste(gene, collapse = ', '))
  }

  # Generate plot for Z
  p1 <- plot_Z(tree)

  # Generate plot for Ps
  p2 <- ggplotify::as.grob(function(){
    graphics::par(mar = c(1, 1, 0, 5))
    P <- tree$Ps
    graphics::image(1:nrow(P), 1:N, axes = FALSE, ylab = "", xlab = "",
                    as.matrix(P[,1:N]), breaks = 0:100/100,
                    col = viridis::turbo(100,begin=0.05,end=0.95),useRaster=T)
    graphics::axis(4, at = 1:N, colnames(P)[1:N], cex.axis = 1, las = 1, tick = FALSE,
                   line=-1)
    graphics::abline(#h = seq(0.5, N + 0.5, 1),
      v = seq(0.5, nrow(P) + 0.5,
              1), col = "grey")
    if(N<=15){
      for (i in 1:nrow(P)) {
        for (j in 1:N) {
          txt.temp <- sprintf("%0.3f", P[i, j])
          graphics::text(i, j, txt.temp, cex = 0.7, col = "white")
        }
      }
    }
  })

  # Generate plot for Pb
  p3 <- ggplotify::as.grob(function(){
    #bottom, left, top and right margins
    graphics::par(mar = c(1, 1, 1, 5))
    P <- tree$Pb
    graphics::image(x=1:nrow(P), y=1:S, axes = FALSE, ylab = "", xlab = "",
                    z=as.matrix(P[,1:S]), breaks = 0:100/100,
                    col = viridis::turbo(100,begin=0.05,end=0.95))
    graphics::axis(4, at = 1:S, colnames(P)[1:S], cex.axis = 1, las = 1, tick = FALSE,
                   line=-1)
    graphics::abline(h = seq(0.5, S + 0.5, 1), v = seq(0.5, nrow(P) + 0.5,
                                                       1), col = "grey")
    for (i in 1:nrow(P)) {
      for (j in 1:S) {
        txt.temp <- sprintf("%0.3f",P[i, j])
        if (P[i, j] <= 0.10 | P[i, j] >= 0.90) {
          graphics::text(i, j, txt.temp, cex = 1, col = "white")
        } else {
          graphics::text(i, j, txt.temp, cex = 1)
        }
      }
    }
  })


  ggp<-ggplot2::ggplot(data.frame(x = seq(0.1,1,0.1), y = seq(0.1,1,0.1), z=seq(0.1,1,0.1)),
              ggplot2::aes(x = x, y = y, fill=z)) +
    ggplot2::geom_tile() + viridis::scale_fill_viridis(option="turbo", limits=c(0,1)) +
    ggplot2::theme(legend.position = c(0.45,.5),
          legend.direction="horizontal",
          legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(1.5, "cm"))+labs(fill = "Prob")

  ggp_legend <- cowplot::get_legend(ggp) # Save legend
  grid::grid.newpage()  # Draw empty plot window
  grid::grid.draw(ggp_legend) # Draw legend to window
  p4 <- grid::grid.grab(wrap=T,wrap.grobs=T) # Grab legend

  # Put everything together
  lay <- rbind(1,1,1,2,2,2,3,4)
  p <- gridExtra::arrangeGrob(p1,p2,p3,p4,layout_matrix = lay)
  grid::grid.draw(p)

  # if TRUE, save text file containing all mutations along the branches
  if (save.muts){
    utils::write.table(mut.list, file = paste0(outpath,"/", project, "_muts.txt"),
                       col.names = FALSE, row.names = FALSE,
                       quote = FALSE, sep = '\t')
  }
  # if TRUE, save the plot of the best tree
  if (save.plot){
    ggplot2::ggsave(filename=paste0(project,"_plot_tree.jpg"), plot=p, path=outpath,
                    height=7, width=7, units="in", limitsize=F, device="jpg")
  }
  return(mut.list)
}

#-------------------------------------------------------------------------------
# Remove cutoff value for posterior probabilities from Canopy function
#-------------------------------------------------------------------------------
canopy_post_v2 <- function(sampchain, projectname, K, numchain, burnin, thin,
                           optK, C = NULL, post.config.cutoff = NULL)
{
  if (is.null(C)) {
    C <- diag(nrow(sampchain[[1]][[1]][[1]]$cna))
    colnames(C) <- rownames(C) <- rownames(sampchain[[1]][[1]][[1]]$cna)
  }
  sampchaink <- sampchain[[which(K == optK)]]
  numchain <- length(sampchaink)
  samptreenew <- sampchaink[[1]][(burnin + 1):length(sampchaink[[1]])]
  numpostburn <- length(samptreenew)
  temp <- thin * c(1:(numpostburn/thin))
  samptreethin <- samptreenew[temp]
  length(samptreethin)
  for (numi in 2:numchain) {
    samptreenew <- sampchaink[[numi]][(burnin + 1):length(sampchaink[[numi]])]
    numpostburn <- length(samptreenew)
    temp <- thin * c(1:(numpostburn/thin))
    samptreethin <- c(samptreethin, samptreenew[temp])
  }
  samptreethin.lik <- rep(NA, length(samptreethin))
  for (treei in 1:length(samptreethin)) {
    samptreethin.lik[treei] <- samptreethin[[treei]]$likelihood
  }
  samptreethin <- samptreethin[which((rank(-1 * samptreethin.lik,
                                           ties.method = "first")) <= 5 * (length(samptreethin)/numchain))]
  samptreethin.lik <- rep(NA, length(samptreethin))
  for (treei in 1:length(samptreethin)) {
    samptreethin.lik[treei] <- samptreethin[[treei]]$likelihood
  }
  if (!is.null(sampchain[[1]][[1]][[1]]$cna)) {
    for (i in 1:length(samptreethin)) {
      samptreethin[[i]] <- sortcna(samptreethin[[i]], C)
    }
  }
  for (i in 1:length(samptreethin)) {
    samptreethin[[i]]$clonalmut <- getclonalcomposition(samptreethin[[i]])
  }
  config <- rep(NA, length(samptreethin))
  config[1] <- 1
  categ <- 1
  for (i in 2:length(samptreethin)) {
    for (categi in 1:categ) {
      list.a <- samptreethin[[i]]$clonalmut
      list.b <- samptreethin[[which(config == categi)[1]]]$clonalmut
      if ((sum(is.element(list.a, list.b)) == optK) & (sum(is.element(list.b,
                                                                      list.a)) == optK)) {
        config[i] <- categi
      }
    }
    if (is.na(config[i])) {
      config[i] <- categ + 1
      categ <- categ + 1
    }
  }
  z.temp <- (samptreethin.lik - mean(samptreethin.lik))/sd(samptreethin.lik)
  samptreethin <- samptreethin[z.temp <= 1.5 & z.temp >= -1.5]
  samptreethin.lik <- samptreethin.lik[z.temp <= 1.5 & z.temp >=
                                         -1.5]
  config <- config[z.temp <= 1.5 & z.temp >= -1.5]
  config.summary <- matrix(nrow = length(unique(config)), ncol = 3)
  colnames(config.summary) <- c("Configuration", "Post_prob",
                                "Mean_post_lik")
  config.summary[, 1] <- unique(config)
  for (i in 1:nrow(config.summary)) {
    configi <- config.summary[i, 1]
    configi.temp <- which(config == configi)
    config.summary[i, 2] <- round(length(configi.temp)/length(config),
                                  3)
    config.summary[i, 3] <- round(max(samptreethin.lik[which(config ==
                                                               configi)]), 2)
  }

  for (treei in 1:length(samptreethin)) {
    output.tree <- samptreethin[[treei]]
    output.tree.Z <- output.tree$Z[, 2:ncol(output.tree$Z),
                                   drop = FALSE]
    output.tree.P <- apply(output.tree$P[2:nrow(output.tree$P),
                                         , drop = FALSE], 2, function(x) {
                                           x/sum(x)
                                         })
    output.tree$CCF <- output.tree.Z %*% output.tree.P
    samptreethin[[treei]] <- output.tree
  }
  return(list(samptreethin, samptreethin.lik, config, config.summary))
}
