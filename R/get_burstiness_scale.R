#' Estimate Bursting Kinetics with SCALE
#'
#' Use the SCALE methodology \insertCite{Jiang2017}{Canopy2} applied to single-cell
#' gene expression data to estimate bursting kinetics parameters: gene
#' activation rate (\code{alpha}), gene deactivation rate (\code{beta}), and
#' transcription rate (\code{scale}). The function performs library size factor
#' normalization internally. This methodology differs from the BPSC methodology
#' \insertCite{Vu2016}{Canopy2} utilized in \code{get_burstiness_bpsc} in that
#' it estimates the parameters by the method of moments. Moment estimators are
#' useful for large datasets as they do not require convergence. However, they
#' can be problematic when the estimates are negative or do not exist.
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
#'
#' @import Rdpack
#' @export
#'
#' @examples
#' # Load post-processed data for patient GBM10
#' data("GBM10_postproc")
#'
#' # Estimate bursting kinetics for the full dataset
#' param.est<-get_burstiness_scale(counts=GBM10_postproc@featurecounts.qc)
#' param.est
#'
#' @references{
#'   \insertAllCited{}
#' }
get_burstiness_scale<-function(counts){

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

  # Next, convert library sizes into size factors by scaling them so that their mean
  # across cells is unity. This ensures that the normalized values are still on the
  # same scale as the raw counts.
  size.factors<-lib.sizes/mean(lib.sizes)

  # Apply the size factors to the raw counts to produce normalized counts
  normalized.counts<-round(counts*size.factors)

  # Function to compute moment estimators by SCALE
  moment_fun<-function(Y.temp,cellsize=NULL){
    if(is.null(cellsize)){cellsize<-rep(1, length(Y.temp))}
    m1<-sum(Y.temp)/sum(cellsize)
    m2<-sum(Y.temp*(Y.temp-1))/sum(cellsize^2)
    m3<-sum(Y.temp*(Y.temp-1)*(Y.temp-2))/sum(cellsize^3)
    kon.hat<--2*(-m1*m2^2+m1^2*m3)/(-m1*m2^2+2*m1^2*m3-m2*m3)
    koff.hat<-2*(m1^2-m2)*(m1*m2-m3)*(m1*m3-m2^2)/(m1^2*m2-2*m2^2+m1*m3)/(2*m1^2*m3-m1*m2^2-m2*m3)
    s.hat<-(-m1*m2^2+2*m1^2*m3-m2*m3)/(m1^2*m2-2*m2^2+m1*m3)
    kinetic.estimates<-list(kon.hat=round(kon.hat,4),
                            koff.hat=round(koff.hat,4),
                            scale=round(s.hat,2))
    return(kinetic.estimates)
  }

  allelic.kinetics.obj<-apply(normalized.counts, 1, function(x) moment_fun(x))
  alpha <-sapply(1:length(allelic.kinetics.obj),
                 function(x) allelic.kinetics.obj[[x]]$kon.hat)
  beta <-sapply(1:length(allelic.kinetics.obj),
                function(x) allelic.kinetics.obj[[x]]$koff.hat)
  scale <-sapply(1:length(allelic.kinetics.obj),
                function(x) allelic.kinetics.obj[[x]]$scale)

  # Indices of genes that have estimable parameters
  id.g<-which(alpha>0 & beta>0 & scale >0)

  # Percent estimable
  pct.estimable<-length(id.g)/length(alpha)

  # Subset include non-negative values
  alpha = alpha[id.g] #gene activation rate
  beta = beta[id.g] #gene deactivation rate
  scale = scale[id.g] #transcription rate

  return(list(alpha=alpha, beta=beta, scale=scale, id.g=id.g, pct.estimable=pct.estimable))
}
