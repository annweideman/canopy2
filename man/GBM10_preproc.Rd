% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gbm_data.R
\docType{data}
\name{GBM10_preproc}
\alias{GBM10_preproc}
\title{Pre-processed glioblastoma (GBM) data from patient GBM10}
\format{
Labels:
\itemize{
\item Bulk samples (RNA-seq): \if{html}{\out{<br />}}
_EGAR00001456512_GBM10_1.RNA, \if{html}{\out{<br />}}
_EGAR00001456513_GBM10_2.RNA, \if{html}{\out{<br />}}
_EGAR00001456514_GBM10_3.RNA
\item Single-cell samples (RNA-seq): \if{html}{\out{<br />}}
_EGAR00001455959_GBM10_1_SC06.RNA to _EGAR00001456043_GBM10_2_SC94.RNA
\item Blood (normal) samples: \if{html}{\out{<br />}}
_EGAR00001456508_GBM10.DNA_B
\item Bulk samples (DNA WES): \if{html}{\out{<br />}}
_EGAR00001456471_GBM10_1.DNA, \if{html}{\out{<br />}}
_EGAR00001456473_GBM10_2.DNA \if{html}{\out{<br />}}
_EGAR00001456472_GBM10_3.DNA
}

\code{GBM10_preproc}   An object of class \code{S4}. Access internal
objects from slots using @.

\code{GBM10_preproc@alt.qc}   A 47 x 84 matrix of pre-processed alternative
read counts from the bulk and single-cell data for GBM10. The row names denote
the chromosomal position of the point mutations. The column names denote the
bulk RNA-seq, scRNA-seq, bulk DNA WES, and blood samples according to the
above labels from the SNV callset of patient GBM10.

\code{GBM10_preproc@annovar}   A 47 x 39 dataframe containing the
functional annotation of the somatic variants for GBM10. For details regarding
the column headers, we refer readers to the extensive
ANNOVAR \href{https://annovar.openbioinformatics.org/}{documentation} and the
original research article \insertCite{Wang2010}{canopy2}.

\code{GBM10_preproc@featurecounts.qc}   A 63677 x 77 matrix of
pre-processed mapped reads from single-cell gene expression data. The row
names denote the gene names. The column names from the third column onward
denote the scRNA-seq samples from the gene expression data for patient GBM10.

\code{GBM10_preproc@mut}   A 47 x 5 dataframe containing the chromosomal
positions of point mutations for GBM10. The columns consist of: the chromosome
number \code{Chr}, the starting position of the point mutation \code{Start},
the ending position of the point mutation \code{End} (coincides with
\code{Start} because the mutation is a single nucleotide variant), the
reference allele \code{Ref}, and the alternate allele \code{Alt}.

\code{GBM10_preproc@ref.qc}   A 47 x 84 matrix of pre-processed
reference read counts from the bulk and single-cell data for GBM10. The
row names denote the chromosomal position of the point mutations. The
column names denote the bulk RNA-seq, scRNA-seq, bulk DNA WES, and blood
samples from the total read counts of patient GBM10.

\code{GBM10_preproc@summary.qc}   A 14 x 77 matrix of the
\code{featureCounts} output summary that reports assignment of alignments
to genomic features. The row names denote the assignment of the alignments.
The column names denote the scRNA-seq samples from the gene expression data
for patient GBM10.
}
\usage{
GBM10_preproc
}
\description{
Pre-processed single-cell RNA sequencing (scRNA-seq), bulk RNA sequencing
(bulk RNA-seq) and bulk DNA whole exome sequencing data (bulk DNA WES) from
one patient (GBM10) selected from 52 individuals with glioblastoma (GBM).
The authors characterized genomic and expression profiles across 127
multisector or longitudinal specimens with matched tumor/normal samples
\insertCite{Lee2017}{canopy2}. The single-nucleotide variant (SNV) pipeline
utilized to pre-process the data is described in the main text and supplement.
Tumor samples were acquired based on 4 categories: 1) locally adjacent
tumors, 2) multifocal/multicentric tumors, 3) 5-ALA (+/-) tumors and
4) Longitudinal tumors. \if{html}{\out{<br />}}
Accession number: EGAS00001001880
}
\details{
\tabular{lllll}{
   \emph{\strong{Patient}} \tab \emph{\strong{Subtype}} \tab \emph{\strong{# bulk DNA WES}} \tab \emph{\strong{# bulk RNA-seq}} \tab \emph{\strong{# scRNA-seq (pre-processed)}} \cr
   ----- \tab ----- \tab ----- \tab ----- \tab ----- \cr
   GBM10 \tab Blood \tab 1 \tab 0 \tab 0 \cr
    \tab GBM10_1 \tab 1 \tab 1 \tab 31 (29) \cr
    \tab GBM10_2 \tab 1 \tab 1 \tab 54 (48) \cr
    \tab GBM10_3 \tab 1 \tab 1 \tab 0 \cr
    \tab (total) \tab 4 \tab 3 \tab 85 (77) \cr
}
}
\references{
{
\insertAllCited{}
}
}
\keyword{datasets}
