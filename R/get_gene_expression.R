#' Secondary processing of gene expression data
#'
# Further pre-process gene expression data:
# 1. Map Ensembl IDs to HGNC symbols
# 2. Remove any duplicates for Ensembl ID
#'
#' @param featurecounts a gene x sample matrix of pre-processed read counts from
#' gene expression data.
#'
#' @param a positive integer corresponding to the Genome Reference Consortium
#' Human (GRCh) build. This is passed to the biomaRt::useEnsembl() function which can
#' currently only take build 37 (as of June 11, 2023).
#'
#' @return
#' The original matrix with two extra columns (appended as the first two columns)
#' that contain the Ensembl gene IDs and corresponding HGNC symbols.
#'
#' @examples
#' #Load pre-processed data for patient GBM10
#' data("GBM10_preproc")
#' featurecounts.qc0<-GBM10_preproc@featurecounts.qc
#'
#' # Get HGNC symbol associated with ensembl gene id
#' featurecounts.qc<-get_gene_expression(featurecounts=featurecounts.qc0,
#'                                       build=37)
#'
#' head(featurecounts.qc)
#'
#' @export

get_gene_expression <- function(featurecounts, build=37){

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
