#' Breast cancer (BC) Data from Patient BC03
#'
#' QC’d single-cell RNA sequencing (scRNA-seq), bulk RNA sequencing (bulk RNA-seq),
#' and bulk DNA whole exome sequencing data (bulk DNA WES) from one patient
#' (BC03) of 11 individuals (BC01-BC11) with primary breast cancer
#' \insertCite{Chung2017}{Canopy2}. The single-nucleotide variant (SNV) pipeline
#' utilized to QC the data is described in Figure S20 of
#' \insertCite{Wang2021}{Canopy2}. In the main study, the authors sequenced 549
#' primary breast cancer cells and matched bulk tumors and/or pooled cells from
#' 11 patients, including two lymph node metastases (BC03LN, BC07LN). Their
#' final dataset after filtering included 515 single cells and 14 bulk samples.
#' Accession Number: GSE75688
#'
#' @format
#'
#' Labels:
#' * Bulk samples (RNA-seq): SRR2973275 (primary sample), SRR2973276 (LN sample)
#' * Single-cell samples (RNA-seq): SRR2973351 to SRR2973436 and SRR5023442
#' * Blood (normal) samples: SRR3023076
#' * Bulk samples (DNA WES): SRR3023078 (LN sample), SRR3023080 (primary sample)
#'
#' \code{BC03} &nbsp; An object of class \code{S4}. Access internal objects
#' from slots using @.
#'
#' \code{BC03@alt.qc} &nbsp; A matrix with 80 columns and 145 rows. The column
#' names denote the bulk RNA-seq, scRNA-seq, bulk DNA WES, and blood samples
#' according to the above labels from the SNV callset of patient BC03. The row
#' names denote the chromosomal position of the point mutations.
#'
#' \code{BC03@annovar} &nbsp; A dataframe with 39 columns and 145 rows.
#' For details regarding the column headers, we refer readers to the extensive
#' ANNOVAR [documentation](https://annovar.openbioinformatics.org/) and the
#' original research article \insertCite{Wang2010}{Canopy2}.
#'
#' \code{BC03@featurecounts.qc} &nbsp; A matrix with 75 columns and 63677 rows.
#' The column names denote the scRNA-seq samples from the gene expression data
#' for patient BC03. The row names denote the gene names.
#'
#' \code{BC03@mut} &nbsp; A data frame with five variables: the chromosome
#' number \code{Chr}, the starting position of the point mutation \code{Start},
#' the ending position of the point mutation \code{End} (coinicides with
#' \code{Start} because the mutation is a single nucleotide variant), the
#' reference allele \code{Ref}, and the alternate allele \code{Alt}.
#'
#' \code{BC03@ref.qc} &nbsp; A matrix with 80 columns and 145 rows. The
#' column names denote the bulk RNA-seq, scRNA-seq, bulk DNA WES, and blood
#' samples from the total read counts of patient BC03. The row names denote the
#' chromosomal position of the point mutations.
#'
#' \code{BC03@summary.qc} &nbsp; A matrix with 75 columns and 14 rows. The
#' column names denote the scRNA-seq samples from the gene expression data for
#' patient BC03. The row names denote the assignment of the alignments.
#'
#' @details
#'
#' \code{BC03@alt.qc} &nbsp; QC'd alternative read counts from bulk and
#' single-cell data for BC03.
#'
#' \code{BC03@annovar} &nbsp; Functional annotation of somatic variants for BC03.
#'
#' \code{BC03@featurecounts.qc} &nbsp; Mapped reads from single-cell gene
#' expression data for BC03.
#'
#' \code{BC03@mut} &nbsp; Chromosomal positions of point mutations for BC03.
#'
#' \code{BC03@ref.qc} &nbsp; QC'd reference read counts from bulk and single-cell
#' data for BC03.
#'
#' \code{BC03@summary.qc} &nbsp; \code{featureCounts} output summary that reports
#' assignment of alignments to genomic features
#'
#'@details
#'
#'  | ***Patient*** | ***Subtype*** | ***# bulk DNA WES*** | ***# bulk RNA-seq*** | ***# scRNA-seq (QC)***|
#'  |-----	|-----	|-----	|-----	|-----	|
#'  |-----	|-----	|-----	|-----	|-----	|
#'  | BC03  | Blood                     	| 1              	| 0              	| 0                	|
#'  |       | Metastatic lymph node       | 1              	| 1 (pooled bulk RNA) | 55 (45)          	|
#'  |       | Primary tumor               | 1              	| 1 (pooled bulk RNA) | 37 (30)          	|
#'  |       | (total)                   	| 3              	| 2              	| 92 (75)        	|
#'
#' @md
#'
#' @import Rdpack
#' @references{
#'    \insertAllCited{}
#' }
#'
"BC03"

#' Breast cancer (BC) Data from Patient BC07
#'
#' QC’d single-cell RNA sequencing (scRNA-seq), bulk RNA sequencing (bulk RNA-seq),
#' and bulk DNA whole exome sequencing data (bulk DNA WES) from one patient
#' (BC07) of 11 individuals (BC01-BC11) with primary breast cancer
#' \insertCite{Chung2017}{Canopy2}. The single-nucleotide variant (SNV) pipeline
#' utilized to QC the data is described in Figure S20 of
#' \insertCite{Wang2021}{Canopy2}. In the main study, the authors sequenced 549
#' primary breast cancer cells and matched bulk tumors and/or pooled cells from
#' 11 patients, including two lymph node metastases (BC03LN, BC07LN). Their
#' final dataset after filtering included 515 single cells and 14 bulk samples.
#' GEO Accession Number: GSE75688
#'
#' @format
#'
#' Labels:
#' * Bulk samples (RNA-seq): SRR2973277 (primary sample), SRR2973278 (LN sample)
#' * Single-cell samples (RNA-seq): SRR2973437 to SRR2973535 and SRR5023558-SRR5023562
#' * Blood (normal) samples: SRR3023081
#' * Bulk samples (DNA WES): SRR3023082 (LN sample), SRR3023083 (primary sample)
#'
#' \code{BC07} &nbsp; An object of class \code{S4}. Access internal objects
#' from slots using @.
#'
#' \code{BC07@alt.qc} &nbsp; A matrix with 106 columns and 79 rows. The column
#' names denote the bulk RNA-seq, scRNA-seq, bulk DNA WES, and blood samples
#' according to the above labels from the SNV callset of patient BC07. The row
#' names denote the chromosomal position of the point mutations.
#'
#' \code{BC07@annovar} &nbsp; A dataframe with 39 columns and 79 rows.
#' For details regarding the column headers, we refer readers to the extensive
#' ANNOVAR [documentation](https://annovar.openbioinformatics.org/) and the
#' original research article \insertCite{Wang2010}{Canopy2}.
#'
#' \code{BC07@featurecounts.qc} &nbsp; A matrix with 101 columns and 63677 rows.
#' The column names denote the scRNA-seq samples from the gene expression data
#' for patient BC07. The row names denote the gene names.
#'
#' \code{BC07@mut} &nbsp; A data frame with five variables: the chromosome
#' number \code{Chr}, the starting position of the point mutation \code{Start},
#' the ending position of the point mutation \code{End} (coinicides with
#' \code{Start} because the mutation is a single nucleotide variant), the
#' reference allele \code{Ref}, and the alternate allele \code{Alt}.
#'
#' \code{BC07@ref.qc} &nbsp; A matrix with 106 columns and 79 rows. The
#' column names denote the bulk RNA-seq, scRNA-seq, bulk DNA WES, and blood
#' samples from the total read counts of patient BC07. The row names denote the
#' chromosomal position of the point mutations.
#'
#' \code{BC07@summary.qc} &nbsp; A matrix with 101 columns and 14 rows. The
#' column names denote the scRNA-seq samples from the gene expression data for
#' patient BC07. The row names denote the assignment of the alignments.
#'
#' @details
#'
#' \code{BC07@alt.qc} &nbsp; QC'd alternative read counts from bulk and
#' single-cell data for BC07.
#'
#' \code{BC07@annovar} &nbsp; Functional annotation of somatic variants for BC07.
#'
#' \code{BC07@featurecounts.qc} &nbsp; Mapped reads from single-cell gene
#' expression data for BC07.
#'
#' \code{BC07@mut} &nbsp; Chromosomal positions of point mutations for BC07.
#'
#' \code{BC07@ref.qc} &nbsp; QC'd reference read counts from bulk and single-cell
#' data for BC07.
#'
#' \code{BC07@summary.qc} &nbsp; \code{featureCounts} output summary that reports
#' assignment of alignments to genomic features
#'
#'@details
#'
#'  | ***Patient*** | ***Subtype*** | ***# bulk DNA WES*** | ***# bulk RNA-seq*** | ***# scRNA-seq (QC)***|
#'  |-----	|-----	|-----	|-----	|-----	|
#'  |-----	|-----	|-----	|-----	|-----	|
#'  | BC07  | Blood                     	| 1              	| 0              	| 0                	|
#'  |       | Metastatic lymph node       | 1              	| 1 (pooled bulk RNA) | 53 (52)          	|
#'  |       | Primary tumor               | 1              	| 1 (tumor bulk RNA) | 51 (49)          	|
#'  |       | (total)                   	| 3              	| 2              	| 104 (101)        	|
#'
#' @md
#'
#' @import Rdpack
#' @references{
#'    \insertAllCited{}
#' }
#'
"BC07"
