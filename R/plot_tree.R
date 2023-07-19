#' Plot phylogenetic tree
#'
#' Plot tumor phylogeny to include clonal configuration matrix, \code{Z},
#' cell-to-clone assignment matrix, \code{Ps}, and sample-to-clone assignment matrix,
#' \code{Pb}. See \bold{details} for further information.
#'
#' @param tree an object of class \code{pyhlo} representing the phylogenetic
#' tree.
#'
#' @param save.muts a logical indicating whether a text file containing all
#' mutations along the branches should be saved. If \code{TRUE}, must specify
#' \code{project} and \code{outpath}. Defaults to \code{FALSE}.
#'
#' @param save.plot a logical indicating whether a plot of the best tree should be
#' saved. If \code{TRUE}, must specify \code{project} and \code{outpath}.
#' Defaults to \code{FALSE}.
#'
#' @param outpath a string specifying the location at which to save output
#' generated from \code{save.muts} and \code{save.plot}
#'
#' @param project a string specifying the project name to add to the filename
#' for output generated from \code{save.muts} and \code{save.plot}. The final
#' filenames for \code{save.muts} and \code{save.plot} will be "projectname_muts.txt"
#' and "projectname_plot_tree.jpg," respectively..
#'
#' @return
#' Automatically plots the tree to the graphics device and outputs a matrix
#' of clusters of mutations that fall along the tree branches.
#'
#' @details
#' The plot consists of a three-tier structure.
#'
#' The first (top) component is the clonal configuration matrix, \code{Z}.
#' This is used to generate a tree that is arranged upright such that time runs
#' vertically from the root to the tips of the tree. The leftmost branch
#' corresponds to the population of normal cells, while the remaining branches
#' correspond to the mutant cells. Each red text box, other than the root,
#' denotes a cluster of point mutations (single-nucleotide variants).
#' For example, "M1: 2 SNVs" would imply that cluster M1 contains two point
#' mutations, perhaps mutations in genes ABC and XYZ. Finally, the tips of each
#' tree are labeled with the percentage of cells assigned to each subclone.
#'
#' The second (middle) component is the cell-to-clone assignment matrix,
#' \code{Ps}. This is a binary matrix, with row sums of one, where the red ones
#' indicate that the  cell (row) belongs to the clone (column) and the blue
#' zeros indicate otherwise.
#'
#' The third (bottom) component is the sample-to-clone assignment matrix,
#' \code{Pb}. This matrix of proportions, with row sums of one, indicates
#' what fraction of each bulk sample (row) belongs to each clone (column).
#'
#' @examples
#' # Simulate read counts
#' sims.out<-simulate_data(N=15, S=2, M=5, alpha=0.1, beta=1.0, kappa=1, tau=999,
#'                        Ktrue=4, b.mindepth=30, b.maxdepth=50, sc.mindepth=80,
#'                        sc.maxdepth=120, scale=300, seed=8675309)
#'
#' N <- ncol(sims.out$Rs) # Number of cells
#' S <- ncol(sims.out$Rb) # Number of bulk samples
#' M <- nrow(sims.out$Rs) # Number of point mutations
#'
#' # Generate a tree
#' K<-4
#' text = paste0(paste0(paste(paste0("(", 1:(K - 1), ","),
#'              collapse = ""), K), paste(rep(")", (K - 1)), collapse = ""), ";")
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
#' # Plot tree
#' plot_tree(tree)
#'
#' @export

plot_tree <- function(tree, save.muts=F, save.plot=F,
                      outpath=NULL, project=NULL){

  if (!inherits(tree, "phylo")){
    stop("tree must be of class \"phylo\"")
  }
  if (is.null(tree$snv)){
    stop("tree must contain point mutations (SNVs) under slot tree$snv")
  }
  if (is.null(tree$Z)){
    stop("tree must contain clonal configuration matrix, Z, under slot tree$Z")
  }
  if (is.null(tree$Ps)){
    stop("tree must contain cell-to-clone assignment matrix, Ps, under slot tree$Ps")
  }
  if (is.null(tree$Pb)){
    stop("tree must contain sample-to-clone assignment matrix, Pb, under slot tree$Pb")
  }
  if(!is.logical(save.muts)){
    stop("Argument save.muts must be a logical (TRUE or FALSE).")
  }
  if(!is.logical(save.plot)){
    stop("Argument save.plot must be a logical (TRUE or FALSE).")
  }
  if(is.null(project) & (save.muts | save.plot)){
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
  if(!is.null(outpath) & !(save.muts | save.plot)){
    stop("You must set save.muts and/or save.plot to 'TRUE' so
         that the output can be directed to the location specified by \"outpath\".")
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

  snvedge <- rep(NA, nrow(tree$snv))
  for (k in 1:nrow(tree$snv)) {
    snvedge[k] <- intersect(which(tree$edge[, 1] == tree$snv[k, 2]),
                            which(tree$edge[, 2] == tree$snv[k, 3]))
  }

  edge.label <- sort(unique(snvedge))
  mut.branch.mat=matrix(nrow=length(edge.label),ncol=1)
  snv.name=rownames(tree$Z)

  # Generate list of mutation clusters
  for (i in 1:length(edge.label)) {
    gene <- snv.name[which(snvedge == edge.label[i])]
    mut.branch.mat[i,1]=paste0("M", i, ": ", paste(gene, collapse = ', '))
    #mut.list<-c(mut.list,paste(gene,collapse = ', '))
  }

  # Generate plot for Z
  p1 <- plot_Z(tree)

  # Generate plot for Ps
  p2 <- ggplotify::as.grob(function(){
    graphics::par(mar = c(1, 1, 0, 5))
    P <- tree$Ps
    N <- ncol(tree$Ps)
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
    S <- ncol(tree$Pb)
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

  temp.df<-data.frame(x = seq(0.1,1,0.1), y = seq(0.1,1,0.1), z=seq(0.1,1,0.1))
  ggp<-ggplot2::ggplot(temp.df,
                       ggplot2::aes(x = temp.df$x, y = temp.df$y, fill=temp.df$z)) +
    ggplot2::geom_tile() + viridis::scale_fill_viridis(option="turbo", limits=c(0,1)) +
    ggplot2::theme(legend.position = c(0.45,.5),
                   legend.direction="horizontal",
                   legend.key.height = grid::unit(0.5, 'cm'),
                   legend.key.width = grid::unit(1.5, "cm"))+ggplot2::labs(fill = "Prob")

  ggp_legend <- cowplot::get_legend(ggp) # Save legend
  grid::grid.newpage()  # Draw empty plot window
  grid::grid.draw(ggp_legend) # Draw legend to window
  p4 <- grid::grid.grab(wrap=T,wrap.grobs=T) # Grab legend

  # Put everything together
  lay <- rbind(1,1,1,2,2,2,3,4)
  p <- gridExtra::arrangeGrob(p1,p2,p3,p4,layout_matrix = lay)
  grid::grid.draw(p)

  # if TRUE, save the plot of the best tree
  if (save.plot){
    ggplot2::ggsave(filename=paste0(project,"_plot_tree.jpg"), plot=p, path=outpath,
                    height=7, width=7, units="in", limitsize=F, device="jpg")
  }
  return(mut.branch.mat)
}
