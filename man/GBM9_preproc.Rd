% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gbm_data.R
\docType{data}
\name{GBM9_preproc}
\alias{GBM9_preproc}
\title{Pre-processed glioblastoma (GBM) data from patient GBM9}
\format{
Labels:
\itemize{
\item Bulk samples (RNA-seq): \if{html}{\out{<br />}}
_EGAR00001456543_GBM9_1.RNA (initial tumor), \if{html}{\out{<br />}}
_EGAR00001456544_GBM9_2.RNA (initial tumor), \if{html}{\out{<br />}}
_EGAR00001456545_GBM9_R1.RNA (recurrent tumor), \if{html}{\out{<br />}}
_EGAR00001456546_GBM9_R2.RNA (recurrent tumor)
\item Single-cell samples (RNA-seq): \if{html}{\out{<br />}}
_EGAR00001456181_GBM9_1_SC07.RNA to _EGAR00001456269_GBM9_2_SC94.RNA (initial tumor) \if{html}{\out{<br />}}
& \if{html}{\out{<br />}}
_EGAR00001456270_GBM9_R1_SC05.RNA to _EGAR00001456313_GBM9_R1_SC93.RNA (recurrent tumor)
\item Blood (normal) samples: \if{html}{\out{<br />}}
_EGAR00001456508_GBM9.DNA_B
\item Bulk samples (DNA WES): \if{html}{\out{<br />}}
_EGAR00001456467_GBM9_1.DNA (initial tumor), \if{html}{\out{<br />}}
_EGAR00001456468_GBM9_2.DNA (initial tumor), \if{html}{\out{<br />}}
_EGAR00001456469_GBM9_R1.DNA (recurrent tumor), \if{html}{\out{<br />}}
_EGAR00001456470_GBM9_R2.DNA (recurrent tumor)
}

\code{GBM9_preproc}   An object of class \code{S4}. Access internal
objects from slots using @.

\code{GBM9_preproc@alt.qc}   A 143 x 135 matrix of pre-processed alternative
read counts from the bulk and single-cell data for GBM9. The row names denote
the chromosomal position of the point mutations. The column names denote the
bulk RNA-seq, scRNA-seq, bulk DNA WES, and blood samples according to the
above labels from the SNV callset of patient GBM9.

\code{GBM9_preproc@annovar}   A 145 x 39 dataframe containing the
functional annotation of the somatic variants for GBM9. For details regarding
the column headers, we refer readers to the extensive
ANNOVAR \href{https://annovar.openbioinformatics.org/}{documentation} and the
original research article \insertCite{Wang2010}{canopy2}.

\code{GBM9_preproc@featurecounts.qc}   A 63677 x 126 matrix of
pre-processed mapped reads from single-cell gene expression data. The row
names denote the gene names. The column names from the third column onward
denote the scRNA-seq samples from the gene expression data for patient GBM9.

\code{GBM9_preproc@mut}   A 143 x 5 dataframe containing the chromosomal
positions of point mutations for GBM9. The columns consist of: the chromosome
number \code{Chr}, the starting position of the point mutation \code{Start},
the ending position of the point mutation \code{End} (coincides with
\code{Start} because the mutation is a single nucleotide variant), the
reference allele \code{Ref}, and the alternate allele \code{Alt}.

\code{GBM9_preproc@ref.qc}   A 143 x 135 matrix of pre-processed
reference read counts from the bulk and single-cell data for GBM9. The
row names denote the chromosomal position of the point mutations. The
column names denote the bulk RNA-seq, scRNA-seq, bulk DNA WES, and blood
samples from the total read counts of patient GBM9.

\code{GBM9_preproc@summary.qc}   A 14 x 126 matrix of the
\code{featureCounts} output summary that reports assignment of alignments
to genomic features. The row names denote the assignment of the alignments.
The column names denote the scRNA-seq samples from the gene expression data
for patient GBM9.
}
\usage{
GBM9_preproc
}
\description{
Pre-processed single-cell RNA sequencing (scRNA-seq), bulk RNA sequencing
(bulk RNA-seq) and bulk DNA whole exome sequencing data (bulk DNA WES) from
one patient (GBM9) selected from 52 individuals with glioblastoma (GBM).
The authors characterized genomic and expression profiles across 127
multisector or longitudinal specimens with matched tumor/normal samples
\insertCite{Lee2017}{canopy2}. The single-nucleotide variant (SNV) pipeline
utilized to pre-process the data is described in the main text and
supplement. Tumor samples were acquired based on 4 categories: 1) locally
adjacent tumors, 2) multifocal/multicentric tumors, 3) 5-ALA (+/-) tumors and
4) Longitudinal tumors. \if{html}{\out{<br />}}
Accession number: EGAS00001001880
}
\details{
\tabular{lllll}{
   \emph{\strong{Patient}} \tab \emph{\strong{Subtype}} \tab \emph{\strong{# bulk DNA WES}} \tab \emph{\strong{# bulk RNA-seq}} \tab \emph{\strong{# scRNA-seq (pre-processed)}} \cr
   ----- \tab ----- \tab ----- \tab ----- \tab ----- \cr
   GBM9 \tab Blood \tab 1 \tab 0 \tab 0 \cr
    \tab GBM9_1 (initial tumor) \tab 1 \tab 1 \tab 29 (25) \cr
    \tab GBM9_2 (initial tumor) \tab 1 \tab 1 \tab 60 (58) \cr
    \tab GBM9_R1 (recurrent tumor) \tab 1 \tab 1 \tab 44 (43) \cr
    \tab GMB9_R2 (recurrent tumor) \tab 1 \tab 1 \tab 0 \cr
    \tab (total) \tab 5 \tab 4 \tab 133 (126) \cr
}
}
\references{
{
\insertAllCited{}
}
}
\keyword{datasets}
