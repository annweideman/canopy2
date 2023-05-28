#' Estimate Bursting Kinetics with BPSC
#'
#' Use the BPSC methodology \insertCite{Vu2016}{Canopy2} applied to single-cell
#' gene expression data to estimate bursting kinetics parameters: gene
#' activation rate (\code{alpha}), gene deactivation rate (\code{beta}), and transcription
#' rate (\code{scale}). The function performs library size factor normalization internally.
#' This methodology differs from the SCALE methodology \insertCite{Jiang2017}{Canopy2}
#' utilized in \code{get_burstiness_scale} in that it estimates the
#' parameters by MCMC sampling. While BPSC tends to be less computationally efficient
#' than SCALE, we found that it had improved estimability (non-NA or negative
#' values).
#'
#' @param counts a gene x single-cell matrix of pre-processed scRNA-seq read counts.
#' Only include count data; do not include any other features/columns like
#' Ensembl ID or HGNC symbol. There is no need to normalize the counts by size
#' factors, as this is performed by the function before computing the estimates.
#'
#' @return
#' A list containing: \code{alpha}, a vector of length M (mutations) that
#' contains the mutation-specific gene activation rates, \code{beta}, a vector
#' of length M containing the mutation-specific gene deactivation rates,
#' \code{scale}, a vector of length M containing the mutation-specific
#' transcription rates, \code{id.g}, a vector of IDs for which parameters can be
#' estimated, and \code{pct.estimable}, a numeric indicating the percent which
#' could be estimated.
#' @import Rdpack
#' @export
#'
#' @examples
#' # Load post-processed data for patient GBM10
#' data("GBM10_postproc")
#'
#' # Estimate bursting kinetics for the full dataset
#' param.est<-get_burstiness_bpsc(counts=GBM10_postproc@featurecounts.qc)
#' param.est
#'
#' @references{
#'   \insertAllCited{}
#' }
get_burstiness_bpsc<-function(counts){

  # check arguments
  if (!inherits(counts, "matrix")){
    stop("counts must be of class \"matrix\"")
  }
  if(!all(counts==round(counts)) | !is.numeric(counts) |
     !all(counts>=0)){
    stop("counts must contain positive integers")
  }

  # Determine the size factors per gene across cells
  lib.sizes<-rowSums(counts)

  # Convert library sizes into size factors by scaling them so that their mean
  # across cells is unity. This ensures that the normalized values are still on the
  # same scale as the raw counts.
  size.factors<-lib.sizes/mean(lib.sizes)

  # Apply the size factors to the raw counts to produce normalized counts
  normalized.counts<-counts*size.factors

  # Estimate parameters from normalized data
  mat.BP<-BPSC::estimateBPMatrix(normalized.counts, para.num=4, estIntPar=F)

  # Parameters (estimated per gene) are every 15th item in the list, e.g., 1, 16, 31, ...
  theta.g<-lapply(seq(1,length(mat.BP$bp.model.list),15),
                  function(x) mat.BP$bp.model.list[[x]])

  # Indices of genes (not all parameters can be estimated, e.g., all zeros or mostly zeros)
  id.g<-sapply(seq(15,length(mat.BP$bp.model.list),15),
               function(x) mat.BP$bp.model.list[[x]])

  # Percent estimable
  pct.estimable<-length(id.g)/nrow(normalized.counts)

  # alpha - gene activation rate
  alpha<-sapply(1:length(theta.g), function(x) theta.g[[x]][1])
  # beta - gene deactivation rate
  beta<-sapply(1:length(theta.g), function(x) theta.g[[x]][2])
  # scale - transcription rate
  scale<-sapply(1:length(theta.g), function(x) theta.g[[x]][3])

  return(list(alpha=alpha, beta=beta, scale=scale,
              id.g=id.g, pct.estimable=pct.estimable))
}
