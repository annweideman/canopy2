#' Map chromosome positions to HGNC symbols
#'
#' Use biomaRt to map chromosome positions to HGNC symbols and return original
#' read counts with appended HGNC symbols.
#'
#' @param reads a matrix with row names of the form chr(num):(pos) where num is a
#' positive integer representing the chromosome number and pos is the positive
#' integer representing the variant position, e.g., chr1:75190429, and
#' column names representing the samples (i.e., single-cell or bulk).
#'
#' @return
#' The original matrix with an extra column (appended as the first column) that
#' contains the HGNC symbols corresponding to the chromosome positions.
#'
#' @examples
#' #Load post-processed data for patient GBM10
#' data("GBM10_postproc")
#'
#' # Map chromosome positions to gene IDs
#' Rs<-getBM_fun(Rs)
#'
#' head(Rs)
#'
#' @export

get_BM_fun<-function(reads){

ensembl<-biomaRt::useEnsembl(biomart="ensembl",
                             dataset="hsapiens_gene_ensembl",
                             GRCh=37)
# Chromosome number
chr.num<-gsub("^[^-]+r|:.*", "", rownames(reads))
# Chromosome start position
start.pos<-end.pos<-gsub(".*:", "", rownames(reads))

temp.df<-data.frame()
for(i in 1:length(chr.num)){
  # Get gene names associated with position of mutation
  # Will output a range an upstream/downstream genes, choose very first range
  getBM.out<-biomaRt::getBM(c("ensembl_gene_id","hgnc_symbol","chromosome_name",
                              "start_position","end_position"),
                            filters = c("chromosome_name", "start"),
                            values = list(chr.num[i],start.pos[i]), mart = ensembl)[1,]
  temp.df<-rbind(temp.df,getBM.out)
}

# Append gene symbol to read counts as first column
reads<-cbind.data.frame(temp.df$hgnc_symbol,reads)
colnames(reads)[1]<-"hgnc_symbol"
reads
}
