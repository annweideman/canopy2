% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gbm_data.R
\docType{data}
\name{GBM2_postproc}
\alias{GBM2_postproc}
\title{Post-processed breast cancer (BC) data from patient GBM2}
\format{
Labels:
\itemize{
\item Single-cell samples (RNA-seq): \if{html}{\out{<br />}}
_EGAR00001456045_GBM2_1_SC02.RNA to _EGAR00001456180_GBM2_2_SC95.RNA
\item Blood (normal) samples: \if{html}{\out{<br />}}
_EGAR00001456500_GBM2.DNA_B
\item Bulk samples (DNA WES): \if{html}{\out{<br />}}
_EGAR00001456446_GBM2_1.DNA, \if{html}{\out{<br />}}
_EGAR00001456447_GBM2_2.DNA
}

\code{GBM2_postproc}   An object of class \code{S4}. Access internal
objects from slots using @.

\code{GBM2_postproc@annovar.qc}   A 317 x 39 dataframe containing the
post-processed functional annotations of the somatic variants. For details
regarding the column headers, we refer readers to the extensive
ANNOVAR \href{https://annovar.openbioinformatics.org/}{documentation} and the
original research article \insertCite{Wang2010}{canopy2}. Only those
annotations corresponding to mutations with at least one non-zero count
across all cells were retained.

\code{GBM2_postproc@featurecounts.qc}   A 63677 x 123 matrix of
post-processed mapped reads from the single-cell gene expression data.
The row names denote the gene names. The column names from the third column
onward denote the scRNA-seq samples from the gene expression data for patient
GBM2. Only those cells that have at least one non-zero count across genes
were retained.

\code{GBM2_postproc@param.est}   A list of 5 elements: \code{alpha},
\code{beta}, \code{scale}, \code{id.g}, and \code{pct.estimable}. The first
three are of length equal to the number of mutations where \code{alpha} denotes the
mutation-specific gene activation rates, \code{beta} denotes the mutation-specific
gene deactivation rates, and \code{scale} denotes the mutation-specific
transcription rates. The ids in \code{id.g} denote the positions corresponding
to the mutations that were estimable (returning non-NA or non-negative values),
and \code{pct.estimable} denotes the percentage that were estimable.

\code{GBM2_postproc@Rb}   A 296 x 2 matrix of post-processed alternative
read counts from the bulk DNA WES data. The row names denote the chromosomal
position of the point mutations. The column names denote the bulk samples
according to the above labels from the SNV callset of patient GBM2.
Post-processing involved removing mutations (rows) with zero alternative
read counts across all cells from the single-cell data.

\code{GBM2_postproc@Rs}   A 296 x 121 matrix of post-processed alternative
read counts from the scRNA-seq data. The row names denote the chromosomal
position of the point mutations. The column names denote the scRNA-seq samples
according to the above labels from the SNV callset of patient GBM2.
Post-processing involved i) removing mutations (rows) with zero alternative
read counts across all cells and ii) removing cells (columns) with zero total
(benign + mutated) read counts across all mutations.

\code{GBM2_postproc@Xb}   A 296 x 2 matrix of post-processed total read
counts (benign + mutated) from the bulk DNA WES data. The row names denote
the chromosomal position of the point mutations. The column names denote the
bulk samples according to the above labels from the SNV callset of patient
GBM2. Post-processing involved removing mutations (rows) with zero
alternative read counts across all cells from the single-cell data.

\code{GBM2_postproc@Xs}   A 296 x 121 matrix of post-processed total read
counts (benign + mutated) from the scRNA-seq data. The row names denote the
chromosomal position of the point mutations. The column names denote the
scRNA-seq samples according to the above labels from the SNV callset of
patient GBM2. Post-processing involved i) removing mutations (rows) with zero
alternative read counts across all cells and ii) removing cells (columns)
with zero total (benign + mutated) read counts across all mutations.
}
\usage{
GBM2_postproc
}
\description{
Some further data munging of GBM2_preproc to get GBM2_postproc, which contains
the finalized read counts for Canopy2 functions and vignette. \if{html}{\out{<br />}}
Accession number: EGAS00001001880
}
\details{
\tabular{lllll}{
   \emph{\strong{Patient}} \tab \emph{\strong{Subtype}} \tab \emph{\strong{# bulk DNA WES}} \tab \emph{\strong{# bulk RNA-seq}} \tab \emph{\strong{# scRNA-seq (post-processed)}} \cr
   ----- \tab ----- \tab ----- \tab ----- \tab ----- \cr
   GBM2 \tab Blood \tab 1 \tab 0 \tab 0 \cr
    \tab GBM2_1 \tab 1 \tab 0 \tab 74 (67) \cr
    \tab GBM2_2 \tab 1 \tab 0 \tab 63 (54) \cr
    \tab (total) \tab 3 \tab 0 \tab 137 (121) \cr
}
}
\references{
{
\insertAllCited{}
}
}
\keyword{datasets}
