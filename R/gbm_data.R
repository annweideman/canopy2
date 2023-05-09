#' Glioblastoma (GBM) Data from Patient GBM2
#'
#' QC'd single-cell RNA sequencing (scRNA-seq) and bulk DNA whole exome
#' sequencing data (bulk DNA WES) from one patient (GBM2) selected from
#' 52 individuals with glioblastoma (GBM). The authors characterized genomic and
#' expression profiles across 127 multisector or longitudinal specimens with
#' matched tumor/normal samples \insertCite{Lee2017}{Canopy2}. The
#' single-nucleotide variant (SNV) pipeline utilized to QC the data is described
#' in Figure S20 of \insertCite{Wang2021}{Canopy2}. Tumor samples
#' were acquired based on 4 categories: 1) locally adjacent tumors,
#' 2) multifocal/multicentric tumors, 3) 5-ALA (+/-) tumors and
#' 4) Longitudinal tumors. Accession number: EGAS00001001880
#'
#' @format
#'
#' Labels:
#' * Single-cell samples (RNA-seq): _EGAR00001456045_GBM2_1_SC02.RNA to
#'                                  _EGAR00001456180_GBM2_2_SC95.RNA
#' * Blood (normal) samples: _EGAR00001456500_GBM2.DNA_B
#' * Bulk samples (DNA WES): _EGAR00001456446_GBM2_1.DNA, _EGAR00001456447_GBM2_2.DNA
#'
#' \code{GBM2} &nbsp; An object of class \code{S4}. Access internal objects
#' from slots using @.
#'
#' \code{GBM2@alt.qc} &nbsp; A matrix with 124 columns and 653 rows. The column
#' names denote the scRNA-seq samples (\code{.RNA}), bulk DNA WES samples
#' (\code{.DNA}), and blood samples (\code{.DNA_B}) from the SNV callset of
#' patient GBM2. The row names denote the chromosomal position of the point
#' mutations.
#'
#' \code{GBM2@annovar} &nbsp; A dataframe with 39 columns and 653 rows.
#' For details regarding the column headers, we refer readers to the extensive
#' [ANNOVAR](https://annovar.openbioinformatics.org/) documentation and the
#' original research article \insertCite{Wang2010}{Canopy2}.
#'
#' \code{GBM2@featurecounts.qc} &nbsp; A matrix with 121 columns and 63677 rows.
#' The column names denote the scRNA-seq samples from the gene expression data
#' for patient GBM2. The row names denote the gene names.
#'
#' \code{GBM2@mut} &nbsp; A data frame with five variables: the chromosome
#' number \code{Chr}, the starting position of the point mutation \code{Start},
#' the ending position of the point mutation \code{End} (coincideswith
#' \code{Start} because the mutation is a single nucleotide variant), the
#' reference allele \code{Ref}, and the alternate allele \code{Alt}.
#'
#' \code{GBM2@ref.qc} &nbsp; A matrix with 124 columns and 653 rows. The
#' column names denote the scRNA-seq samples (\code{.RNA}), bulk DNA WES samples
#' (\code{.DNA}), and blood samples (\code{.DNA_B}) from the total read counts
#' of patient GBM2. The row names denote the chromosomal position of the point
#' mutations.
#'
#' \code{GBM2@summary.qc} &nbsp; A matrix with 121 columns and 14 rows. The
#' column names denote the scRNA-seq samples from the gene expression data for
#' patient GBM2. The row names denote the assignment of the alignments.
#'
#' @details
#'
#' \code{GBM2@alt.qc} &nbsp; Preprocessed alternative read counts from bulk DNA
#' WES and scRNA-seq for GMB2.
#'
#' \code{GBM2@annovar} &nbsp; Functional annotation of somatic variants for GBM2.
#'
#' \code{GBM2@featurecounts.qc} &nbsp; Mapped reads from single-cell gene
#' expression data for GBM2.
#'
#' \code{GBM2@mut} &nbsp; Chromosomal positions of point mutations for GBM2.
#'
#' \code{GBM2@ref.qc} &nbsp; Preprocessed reference read counts from bulk DNA WES and
#' single-cell RNA-seq for GBM2.
#'
#' \code{GBM2@summary.qc} &nbsp; \code{featureCounts} output summary that reports
#' assignment of alignments to genomic features
#'
#'@details
#'
#'  | ***Patient*** | ***Subtype*** | ***# bulk DNA WES*** | ***# bulk RNA-seq*** | ***# scRNA-seq (QC)***|
#'  |-----	|-----	|-----	|-----	|-----	|
#'  |-----	|-----	|-----	|-----	|-----	|
#'  | GBM2  | Blood                     	| 1              	| 0              	| 0                	|
#'  |       | GBM2_1                    	| 1              	| 0              	| 74 (67)          	|
#'  |       | GBM2_2                    	| 1              	| 0              	| 63 (54)          	|
#'  |       | (total)                   	| 3              	| 0              	| 137 (121)        	|
#'
#' @md
#'
#' @import Rdpack
#' @references{
#'    \insertAllCited{}
#' }
#'
"GBM2"

#' Glioblastoma (GBM) Data from Patient GBM9
#'
#' QC'd single-cell RNA sequencing (scRNA-seq), bulk RNA sequencing (bulk RNA-seq)
#' and bulk DNA whole exome sequencing data (bulk DNA WES) from one patient (GBM9)
#' selected from 52 individuals with glioblastoma (GBM). The authors characterized
#' genomic and expression profiles across 127 multisector or longitudinal specimens
#' with matched tumor/normal samples \insertCite{Lee2017}{Canopy2}. The
#' single-nucleotide variant (SNV) pipeline utilized to QC the data is described
#' in Figure S20 of \insertCite{Wang2021}{Canopy2}. Tumor samples
#' were acquired based on 4 categories: 1) locally adjacent tumors,
#' 2) multifocal/multicentric tumors, 3) 5-ALA (+/-) tumors and
#' 4) Longitudinal tumors. Accession number: EGAS00001001880
#'
#' @format
#'
#' Labels:
#' * Bulk samples (RNA-seq): _EGAR00001456543_GBM9_1.RNA (initial tumor),
#'                           _EGAR00001456544_GBM9_2.RNA (initial tumor),
#'                           _EGAR00001456545_GBM9_R1.RNA (recurrent tumor),
#'                           _EGAR00001456546_GBM9_R2.RNA (recurrent tumor)
#' * Single-cell samples (RNA-seq): _EGAR00001456181_GBM9_1_SC07.RNA to
#'                                  _EGAR00001456269_GBM9_2_SC94.RNA (initial tumor)
#'                                  &
#'                                  _EGAR00001456270_GBM9_R1_SC05.RNA to
#'                                  _EGAR00001456313_GBM9_R1_SC93.RNA (recurrent tumor)
#' * Blood (normal) samples: _EGAR00001456508_GBM9.DNA_B
#' * Bulk samples (DNA WES): _EGAR00001456467_GBM9_1.DNA (initial tumor),
#'                           _EGAR00001456468_GBM9_2.DNA (initial tumor),
#'                           _EGAR00001456469_GBM9_R1.DNA (recurrent tumor),
#'                           _EGAR00001456470_GBM9_R2.DNA (recurrent tumor)
#'
#' \code{GBM9} &nbsp; An object of class \code{S4}. Access internal objects
#' from slots using @.
#'
#' \code{GBM9@alt.qc} &nbsp; A matrix with 135 columns and 143 rows. The column
#' names denote the scRNA-seq and bulk RNA-seq samples (\code{.RNA}),
#' bulk DNA WES samples (\code{.DNA}), and blood samples (\code{.DNA_B}) from
#' the SNV callset of patient GBM9. The row names denote the chromosomal
#' position of the point mutations.
#'
#' \code{GBM9@annovar} &nbsp; A dataframe with 39 columns and 143 rows.
#' For details regarding the column headers, we refer readers to the extensive
#' [ANNOVAR](https://annovar.openbioinformatics.org/) documentation and the
#' original research article \insertCite{Wang2010}{Canopy2}.
#'
#' \code{GBM9@featurecounts.qc} &nbsp; A matrix with 126 columns and 63677 rows.
#' The column names denote the scRNA-seq samples from the gene expression data
#' for patient GBM9. The row names denote the gene names.
#'
#' \code{GBM9@mut} &nbsp; A data frame with five variables: the chromosome
#' number \code{Chr}, the starting position of the point mutation \code{Start},
#' the ending position of the point mutation \code{End} (coincideswith
#' \code{Start} because the mutation is a single nucleotide variant), the
#' reference allele \code{Ref}, and the alternate allele \code{Alt}.
#'
#' \code{GBM9@ref.qc} &nbsp; A matrix with 135 columns and 143 rows. The
#' column names denote the scRNA-seq and bulk RNA-seq samples (\code{.RNA}),
#' bulk DNA WES samples (\code{.DNA}), and blood samples (\code{.DNA_B})
#' from the total read counts of patient GBM9. The row names denote the
#' chromosomal position of the point mutations.
#'
#' \code{GBM9@summary.qc} &nbsp; A matrix with 126 columns and 14 rows. The
#' column names denote the scRNA-seq samples from the gene expression data for
#' patient GBM9. The row names denote the assignment of the alignments.
#'
#' @details
#'
#' \code{GBM9@alt.qc} &nbsp; Preprocessed alternative read counts from bulk DNA
#' WES and scRNA-seq for GMB9.
#'
#' \code{GBM9@annovar} &nbsp; Functional annotation of somatic variants for GBM9.
#'
#' \code{GBM9@featurecounts.qc} &nbsp; Mapped reads from single-cell gene
#' expression data for GBM9.
#'
#' \code{GBM9@mut} &nbsp; Chromosomal positions of point mutations for GBM9.
#'
#' \code{GBM9@ref.qc} &nbsp; Preprocessed reference read counts from bulk DNA WES and
#' single-cell RNA-seq for GBM9.
#'
#' \code{GBM9@summary.qc} &nbsp; \code{featureCounts} output summary that reports
#' assignment of alignments to genomic features
#'
#'@details
#'
#'  | ***Patient*** | ***Subtype*** | ***# bulk DNA WES*** | ***# bulk RNA-seq*** | ***# scRNA-seq (QC)***|
#'  |-----	|-----	|-----	|-----	|-----	|
#'  |-----	|-----	|-----	|-----	|-----	|
#'  | GBM9	| Blood                     	| 1              	| 0              	| 0                	|
#'  |       | GBM9_1 (initial tumor)    	| 1              	| 1              	| 29 (25)          	|
#'  |       | GBM9_2 (initial tumor)    	| 1              	| 1              	| 60 (58)          	|
#'  |       | GBM9_R1 (recurrent tumor) 	| 1              	| 1              	| 44 (43)          	|
#'  |       | GMB9_R2 (recurrent tumor) 	| 1              	| 1              	| 0                	|
#'  |       | (total)                   	| 5              	| 4              	| 133 (126)        	|
#'
#' @md
#'
#' @import Rdpack
#' @references{
#'    \insertAllCited{}
#' }
#'
"GBM9"

#' Glioblastoma (GBM) Data from Patient GBM10
#'
#' QC'd single-cell RNA sequencing (scRNA-seq), bulk RNA sequencing (bulk RNA-seq)
#' and bulk DNA whole exome sequencing data (bulk DNA WES) from one patient (GBM10)
#' selected from 52 individuals with glioblastoma (GBM). The authors characterized
#' genomic and expression profiles across 127 multisector or longitudinal specimens
#' with matched tumor/normal samples \insertCite{Lee2017}{Canopy2}. The
#' single-nucleotide variant (SNV) pipeline utilized to QC the data is described
#' in Figure S20 of \insertCite{Wang2021}{Canopy2}. Tumor samples were acquired
#' based on 4 categories: 1) locally adjacent tumors, 2) multifocal/multicentric
#' tumors, 3) 5-ALA (+/-) tumors and 4) Longitudinal tumors.
#' Accession number: EGAS00001001880
#'
#' @format
#'
#' Labels:
#' * Bulk samples (RNA-seq): _EGAR00001456512_GBM10_1.RNA,
#'                           _EGAR00001456513_GBM10_2.RNA,
#'                           _EGAR00001456514_GBM10_3.RNA
#' * Single-cell samples (RNA-seq): _EGAR00001455959_GBM10_1_SC06.RNA to
#'                                  _EGAR00001456043_GBM10_2_SC94.RNA
#' * Blood (normal) samples: _EGAR00001456508_GBM9.DNA_B
#' * Bulk samples (DNA WES): _EGAR00001456471_GBM10_1.DNA,
#'                           _EGAR00001456473_GBM10_2.DNA
#'                           _EGAR00001456472_GBM10_3.DNA
#'
#' @format
#'
#' \code{GBM10} &nbsp; An object of class \code{S4}. Access internal objects
#' from slots using @.
#'
#' \code{GBM10@alt.qc} &nbsp; A matrix with 84 columns and 47 rows. The column
#' names denote the scRNA-seq and bulk RNA-seq samples (\code{.RNA}),
#' bulk DNA WES samples (\code{.DNA}), and blood samples (\code{.DNA_B})
#' from the SNV callset of patient GBM10. The row names denote the chromosomal
#' position of the point mutations.
#'
#' \code{GBM10@annovar} &nbsp; A dataframe with 39 columns and 47 rows.
#' For details regarding the column headers, we refer readers to the extensive
#' [ANNOVAR](https://annovar.openbioinformatics.org/) documentation and the
#' original research article \insertCite{Wang2010}{Canopy2}.
#'
#' \code{GBM10@featurecounts.qc} &nbsp; A matrix with 77 columns and 63677 rows.
#' The column names denote the scRNA-seq samples from the gene expression data
#' for patient GBM10. The row names denote the gene names.
#'
#' \code{GBM10@mut} &nbsp; A data frame with five variables: the chromosome
#' number \code{Chr}, the starting position of the point mutation \code{Start},
#' the ending position of the point mutation \code{End} (coincideswith
#' \code{Start} because the mutation is a single nucleotide variant), the
#' reference allele \code{Ref}, and the alternate allele \code{Alt}.
#'
#' \code{GBM10@ref.qc} &nbsp; A matrix with 84 columns and 47 rows. The
#' column names denote the scRNA-seq and bulk RNA-seq samples (\code{.RNA}),
#' bulk DNA WES samples (\code{.DNA}), and blood samples (\code{.DNA_B}) from
#' the total read counts of patient GBM10. The row names denote the chromosomal
#' position of the point mutations.
#'
#' \code{GBM10@summary.qc} &nbsp; A matrix with 77 columns and 14 rows. The
#' column names denote the scRNA-seq samples from the gene expression data for
#' patient GBM10. The row names denote the assignment of the alignments.
#'
#' @details
#'
#' \code{GBM10@alt.qc} &nbsp; Preprocessed alternative read counts from bulk DNA
#' WES and scRNA-seq for GMB10.
#'
#' \code{GBM10@annovar} &nbsp; Functional annotation of somatic variants for GBM10.
#'
#' \code{GBM10@featurecounts.qc} &nbsp; Mapped reads from single-cell gene
#' expression data for GBM10.
#'
#' \code{GBM10@mut} &nbsp; Chromosomal positions of point mutations for GBM10.
#'
#' \code{GBM10@ref.qc} &nbsp; Preprocessed reference read counts from bulk DNA WES and
#' single-cell RNA-seq for GBM10.
#'
#' \code{GBM10@summary.qc} &nbsp; \code{featureCounts} output summary that reports
#' assignment of alignments to genomic features
#'
#'@details
#'
#'  | ***Patient*** | ***Subtype*** | ***# bulk DNA WES*** | ***# bulk RNA-seq*** | ***# scRNA-seq (QC)***|
#'  |-----	|-----	|-----	|-----	|-----	|
#'  |-----	|-----	|-----	|-----	|-----	|
#'  | GBM10 | Blood                     	| 1              	| 0              	| 0                	|
#'  |       | GBM10_1                   	| 1              	| 1              	| 31 (29)          	|
#'  |       | GBM10_2                   	| 1              	| 1              	| 54 (48)          	|
#'  |       | GBM10_3                   	| 1              	| 1              	| 0                	|
#'  |       | (total)                   	| 4              	| 3              	| 85 (77)          	|
#'
#' @md
#'
#' @import Rdpack
#' @references{
#'    \insertAllCited{}
#' }
#'
"GBM10"
