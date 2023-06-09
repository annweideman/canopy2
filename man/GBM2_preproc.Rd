% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gbm_data.R
\docType{data}
\name{GBM2_preproc}
\alias{GBM2_preproc}
\title{Pre-processed glioblastoma (GBM) Data from patient GBM2}
\format{
Labels:
\itemize{
\item Single-cell samples (RNA-seq): \if{html}{\out{<br />}}
_EGAR00001456045_GBM2_1_SC02.RNA to _EGAR00001456180_GBM2_2_SC95.RNA
\item Blood (normal) samples: \if{html}{\out{<br />}}
_EGAR00001456500_GBM2.DNA_B \if{html}{\out{<br />}}
\item Bulk samples (DNA WES): \if{html}{\out{<br />}}
_EGAR00001456446_GBM2_1.DNA,\if{html}{\out{<br />}}
_EGAR00001456447_GBM2_2.DNA
}

\code{GBM2_preproc}   An object of class \code{S4}. Access internal
objects from slots using @.

\code{GBM2_preproc@alt.qc}   A 653 x 124 matrix of pre-processed
alternative read counts from the bulk and single-cell data for GBM2.
The row names denote the chromosomal position of the point mutations.
The column names denote the scRNA-seq, bulk DNA WES, and blood samples
according to the above labels from the SNV callset of patient GBM2.

\code{GBM2_preproc@annovar}   A 653 x 39 dataframe containing the
functional annotation of the somatic variants for GBM2. For details regarding
the column headers, we refer readers to the extensive
ANNOVAR \href{https://annovar.openbioinformatics.org/}{documentation} and the
original research article \insertCite{Wang2010}{canopy2}.

\code{GBM2_preproc@featurecounts.qc}   A 63677 x 121 matrix of
pre-processed mapped reads from single-cell gene expression data. The row
names denote the gene names. The column names from the third column onward
denote the scRNA-seq samples from the gene expression data for patient GBM2.

\code{GBM2_preproc@mut}   A 653 x 5 dataframe containing the chromosomal
positions of point mutations for GBM2. The columns consist of: the chromosome
number \code{Chr}, the starting position of the point mutation \code{Start},
the ending position of the point mutation \code{End} (coincides with
\code{Start} because the mutation is a single nucleotide variant), the
reference allele \code{Ref}, and the alternate allele \code{Alt}.

\code{GBM2_preproc@ref.qc}   A 653 x 124 matrix of pre-processed
reference read counts from the bulk and single-cell data for GBM2. The
row names denote the chromosomal position of the point mutations. The
column names denote the scRNA-seq, bulk DNA WES, and blood samples from the
total read counts of patient GBM2.

\code{GBM2_preproc@summary.qc}   A 14 x 121 matrix of the
\code{featureCounts} output summary that reports assignment of alignments
to genomic features. The row names denote the assignment of the alignments.
The column names denote the scRNA-seq samples from the gene expression data
for patient GBM2.
}
\usage{
GBM2_preproc
}
\description{
Pre-processed single-cell RNA sequencing (scRNA-seq) and bulk DNA whole exome
sequencing data (bulk DNA WES) from one patient (GBM2) selected from
52 individuals with glioblastoma (GBM). The authors characterized genomic and
expression profiles across 127 multisector or longitudinal specimens with
matched tumor/normal samples \insertCite{Lee2017}{canopy2}. The
single-nucleotide variant (SNV) pipeline utilized to pre-process the data is
described in the main text and supplement. Tumor samples were acquired based
on 4 categories: 1) locally adjacent tumors, 2) multifocal/multicentric
tumors, 3) 5-ALA (+/-) tumors and 4) Longitudinal tumors. \if{html}{\out{<br />}}
Accession number: EGAS00001001880
}
\details{
\tabular{lllll}{
   \emph{\strong{Patient}} \tab \emph{\strong{Subtype}} \tab \emph{\strong{# bulk DNA WES}} \tab \emph{\strong{# bulk RNA-seq}} \tab \emph{\strong{# scRNA-seq (pre-processed)}} \cr
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
