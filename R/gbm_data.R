#' Pre-processed glioblastoma (GBM) Data from patient GBM2
#'
#' Pre-processed single-cell RNA sequencing (scRNA-seq) and bulk DNA whole exome
#' sequencing data (bulk DNA WES) from one patient (GBM2) selected from
#' 52 individuals with glioblastoma (GBM). The authors characterized genomic and
#' expression profiles across 127 multisector or longitudinal specimens with
#' matched tumor/normal samples \insertCite{Lee2017}{Canopy2}. The
#' single-nucleotide variant (SNV) pipeline utilized to pre-process the data is
#' described in the main text and supplement. Tumor samples were acquired based
#' on 4 categories: 1) locally adjacent tumors, 2) multifocal/multicentric
#' tumors, 3) 5-ALA (+/-) tumors and 4) Longitudinal tumors. <br />
#' Accession number: EGAS00001001880
#'
#' @format
#'
#' Labels:
#' * Single-cell samples (RNA-seq): <br />
#' _EGAR00001456045_GBM2_1_SC02.RNA to _EGAR00001456180_GBM2_2_SC95.RNA
#' * Blood (normal) samples: <br />
#' _EGAR00001456500_GBM2.DNA_B <br />
#' * Bulk samples (DNA WES): <br />
#' _EGAR00001456446_GBM2_1.DNA,<br />
#' _EGAR00001456447_GBM2_2.DNA
#'
#' \code{GBM2_preproc} &nbsp; An object of class \code{S4}. Access internal
#' objects from slots using @.
#'
#' \code{GBM2_preproc@alt.qc} &nbsp; A 653 x 124 matrix of pre-processed
#' alternative read counts from the bulk and single-cell data for GBM2.
#' The row names denote the chromosomal position of the point mutations.
#' The column names denote the scRNA-seq, bulk DNA WES, and blood samples
#' according to the above labels from the SNV callset of patient GBM2.
#'
#' \code{GBM2_preproc@annovar} &nbsp; A 653 x 39 dataframe containing the
#' functional annotation of the somatic variants for GBM2. For details regarding
#' the column headers, we refer readers to the extensive
#' ANNOVAR [documentation](https://annovar.openbioinformatics.org/) and the
#' original research article \insertCite{Wang2010}{Canopy2}.
#'
#' \code{GBM2_preproc@featurecounts.qc} &nbsp; A 63677 x 121 matrix of
#' pre-processed mapped reads from single-cell gene expression data. The row
#' names denote the gene names. The column names from the third column onward
#' denote the scRNA-seq samples from the gene expression data for patient GBM2.
#'
#' \code{GBM2_preproc@mut} &nbsp; A 653 x 5 dataframe containing the chromosomal
#' positions of point mutations for GBM2. The columns consist of: the chromosome
#' number \code{Chr}, the starting position of the point mutation \code{Start},
#' the ending position of the point mutation \code{End} (coincides with
#' \code{Start} because the mutation is a single nucleotide variant), the
#' reference allele \code{Ref}, and the alternate allele \code{Alt}.
#'
#' \code{GBM2_preproc@ref.qc} &nbsp; A 653 x 124 matrix of pre-processed
#' reference read counts from the bulk and single-cell data for GBM2. The
#' row names denote the chromosomal position of the point mutations. The
#' column names denote the scRNA-seq, bulk DNA WES, and blood samples from the
#' total read counts of patient GBM2.
#'
#' \code{GBM2_preproc@summary.qc} &nbsp; A 14 x 121 matrix of the
#' \code{featureCounts} output summary that reports assignment of alignments
#' to genomic features. The row names denote the assignment of the alignments.
#' The column names denote the scRNA-seq samples from the gene expression data
#' for patient GBM2.
#'
#'@details
#'
#'  | ***Patient*** | ***Subtype*** | ***# bulk DNA WES*** | ***# bulk RNA-seq*** | ***# scRNA-seq (pre-processed)***|
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
"GBM2_preproc"
#'
#' Post-processed breast cancer (BC) data from patient GBM2
#'
#' Some further data munging of GBM2_preproc to get GBM2_postproc, which contains
#' the finalized read counts for Canopy2 functions and vignette. <br />
#' Accession number: EGAS00001001880
#'
#' @format
#'
#' Labels:
#' * Single-cell samples (RNA-seq): <br />
#' _EGAR00001456045_GBM2_1_SC02.RNA to _EGAR00001456180_GBM2_2_SC95.RNA
#' * Blood (normal) samples: <br />
#'  _EGAR00001456500_GBM2.DNA_B
#' * Bulk samples (DNA WES): <br />
#' _EGAR00001456446_GBM2_1.DNA, <br />
#' _EGAR00001456447_GBM2_2.DNA
#'
#' \code{GBM2_postproc} &nbsp; An object of class \code{S4}. Access internal
#' objects from slots using @.
#'
#' \code{GBM2_postproc@annovar.qc} &nbsp; A 317 x 39 dataframe containing the
#' post-processed functional annotations of the somatic variants. For details
#' regarding the column headers, we refer readers to the extensive
#' ANNOVAR [documentation](https://annovar.openbioinformatics.org/) and the
#' original research article \insertCite{Wang2010}{Canopy2}. Only those
#' annotations corresponding to mutations with at least one non-zero count
#' across all cells were retained.
#'
#' \code{GBM2_postproc@featurecounts.qc} &nbsp; A 63677 x 123 matrix of
#' post-processed mapped reads from the single-cell gene expression data.
#' The row names denote the gene names. The column names from the third column
#' onward denote the scRNA-seq samples from the gene expression data for patient
#' GBM2. Only those cells that have at least one non-zero count across genes
#' were retained.
#'
#' \code{GBM2_postproc@param.est} &nbsp; A list of 5 elements: \code{alpha},
#' \code{beta}, \code{scale}, \code{id.g}, and \code{pct.estimable}. The first
#' three are of length equal to the number of mutations where \code{alpha} denotes the
#' mutation-specific gene activation rates, \code{beta} denotes the mutation-specific
#' gene deactivation rates, and \code{scale} denotes the mutation-specific
#' transcription rates. The ids in \code{id.g} denote the positions corresponding
#' to the mutations that were estimable (returning non-NA or non-negative values),
#' and \code{pct.estimable} denotes the percentage that were estimable.
#'
#' \code{GBM2_postproc@Rb} &nbsp; A 296 x 2 matrix of post-processed alternative
#' read counts from the bulk DNA WES data. The row names denote the chromosomal
#' position of the point mutations. The column names denote the bulk samples
#' according to the above labels from the SNV callset of patient GBM2.
#' Post-processing involved removing mutations (rows) with zero alternative
#' read counts across all cells from the single-cell data.
#'
#' \code{GBM2_postproc@Rs} &nbsp; A 296 x 121 matrix of post-processed alternative
#' read counts from the scRNA-seq data. The row names denote the chromosomal
#' position of the point mutations. The column names denote the scRNA-seq samples
#' according to the above labels from the SNV callset of patient GBM2.
#' Post-processing involved i) removing mutations (rows) with zero alternative
#' read counts across all cells and ii) removing cells (columns) with zero total
#' (benign + mutated) read counts across all mutations.
#'
#' \code{GBM2_postproc@Xb} &nbsp; A 296 x 2 matrix of post-processed total read
#' counts (benign + mutated) from the bulk DNA WES data. The row names denote
#' the chromosomal position of the point mutations. The column names denote the
#' bulk samples according to the above labels from the SNV callset of patient
#' GBM2. Post-processing involved removing mutations (rows) with zero
#' alternative read counts across all cells from the single-cell data.
#'
#' \code{GBM2_postproc@Xs} &nbsp; A 296 x 121 matrix of post-processed total read
#' counts (benign + mutated) from the scRNA-seq data. The row names denote the
#' chromosomal position of the point mutations. The column names denote the
#' scRNA-seq samples according to the above labels from the SNV callset of
#' patient GBM2. Post-processing involved i) removing mutations (rows) with zero
#' alternative read counts across all cells and ii) removing cells (columns)
#' with zero total (benign + mutated) read counts across all mutations.
#'
#'@details
#'
#'  | ***Patient*** | ***Subtype*** | ***# bulk DNA WES*** | ***# bulk RNA-seq*** | ***# scRNA-seq (post-processed)***|
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
"GBM2_postproc"
#'
#' Pre-processed glioblastoma (GBM) data from patient GBM9
#'
#' Pre-processed single-cell RNA sequencing (scRNA-seq), bulk RNA sequencing
#' (bulk RNA-seq) and bulk DNA whole exome sequencing data (bulk DNA WES) from
#' one patient (GBM9) selected from 52 individuals with glioblastoma (GBM).
#' The authors characterized genomic and expression profiles across 127
#' multisector or longitudinal specimens with matched tumor/normal samples
#' \insertCite{Lee2017}{Canopy2}. The single-nucleotide variant (SNV) pipeline
#' utilized to pre-process the data is described in the main text and
#' supplement. Tumor samples were acquired based on 4 categories: 1) locally
#' adjacent tumors, 2) multifocal/multicentric tumors, 3) 5-ALA (+/-) tumors and
#' 4) Longitudinal tumors. <br />
#' Accession number: EGAS00001001880
#'
#' @format
#'
#' Labels:
#' * Bulk samples (RNA-seq): <br />
#' _EGAR00001456543_GBM9_1.RNA (initial tumor), <br />
#' _EGAR00001456544_GBM9_2.RNA (initial tumor), <br />
#' _EGAR00001456545_GBM9_R1.RNA (recurrent tumor), <br />
#' _EGAR00001456546_GBM9_R2.RNA (recurrent tumor)
#' * Single-cell samples (RNA-seq): <br />
#' _EGAR00001456181_GBM9_1_SC07.RNA to _EGAR00001456269_GBM9_2_SC94.RNA (initial tumor) <br />
#' & <br />
#' _EGAR00001456270_GBM9_R1_SC05.RNA to _EGAR00001456313_GBM9_R1_SC93.RNA (recurrent tumor)
#' * Blood (normal) samples: <br />
#' _EGAR00001456508_GBM9.DNA_B
#' * Bulk samples (DNA WES): <br />
#' _EGAR00001456467_GBM9_1.DNA (initial tumor), <br />
#' _EGAR00001456468_GBM9_2.DNA (initial tumor), <br />
#' _EGAR00001456469_GBM9_R1.DNA (recurrent tumor), <br />
#' _EGAR00001456470_GBM9_R2.DNA (recurrent tumor)
#'
#' \code{GBM9_preproc} &nbsp; An object of class \code{S4}. Access internal
#' objects from slots using @.
#'
#' \code{GBM9_preproc@alt.qc} &nbsp; A 143 x 135 matrix of pre-processed alternative
#' read counts from the bulk and single-cell data for GBM9. The row names denote
#' the chromosomal position of the point mutations. The column names denote the
#' bulk RNA-seq, scRNA-seq, bulk DNA WES, and blood samples according to the
#' above labels from the SNV callset of patient GBM9.
#'
#' \code{GBM9_preproc@annovar} &nbsp; A 145 x 39 dataframe containing the
#' functional annotation of the somatic variants for GBM9. For details regarding
#' the column headers, we refer readers to the extensive
#' ANNOVAR [documentation](https://annovar.openbioinformatics.org/) and the
#' original research article \insertCite{Wang2010}{Canopy2}.
#'
#' \code{GBM9_preproc@featurecounts.qc} &nbsp; A 63677 x 126 matrix of
#' pre-processed mapped reads from single-cell gene expression data. The row
#' names denote the gene names. The column names from the third column onward
#' denote the scRNA-seq samples from the gene expression data for patient GBM9.
#'
#' \code{GBM9_preproc@mut} &nbsp; A 143 x 5 dataframe containing the chromosomal
#' positions of point mutations for GBM9. The columns consist of: the chromosome
#' number \code{Chr}, the starting position of the point mutation \code{Start},
#' the ending position of the point mutation \code{End} (coincides with
#' \code{Start} because the mutation is a single nucleotide variant), the
#' reference allele \code{Ref}, and the alternate allele \code{Alt}.
#'
#' \code{GBM9_preproc@ref.qc} &nbsp; A 143 x 135 matrix of pre-processed
#' reference read counts from the bulk and single-cell data for GBM9. The
#' row names denote the chromosomal position of the point mutations. The
#' column names denote the bulk RNA-seq, scRNA-seq, bulk DNA WES, and blood
#' samples from the total read counts of patient GBM9.
#'
#' \code{GBM9_preproc@summary.qc} &nbsp; A 14 x 126 matrix of the
#' \code{featureCounts} output summary that reports assignment of alignments
#' to genomic features. The row names denote the assignment of the alignments.
#' The column names denote the scRNA-seq samples from the gene expression data
#' for patient GBM9.
#'
#'@details
#'
#'  | ***Patient*** | ***Subtype*** | ***# bulk DNA WES*** | ***# bulk RNA-seq*** | ***# scRNA-seq (pre-processed)***|
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
"GBM9_preproc"
#'
#' Post-processed breast cancer (BC) data from patient GBM9
#'
#' Some further data munging of GBM9_preproc to get GBM9_postproc, which contains
#' the finalized read counts for Canopy2 functions and vignette. <br />
#' Accession number: EGAS00001001880
#'
#' @format
#'
#' Labels:
#' * Bulk samples (RNA-seq): <br />
#' _EGAR00001456543_GBM9_1.RNA (initial tumor), <br />
#' _EGAR00001456544_GBM9_2.RNA (initial tumor), <br />
#' _EGAR00001456545_GBM9_R1.RNA (recurrent tumor), <br />
#' _EGAR00001456546_GBM9_R2.RNA (recurrent tumor)
#' * Single-cell samples (RNA-seq): <br />
#' _EGAR00001456181_GBM9_1_SC07.RNA to _EGAR00001456269_GBM9_2_SC94.RNA (initial tumor) <br />
#' & <br />
#' _EGAR00001456270_GBM9_R1_SC05.RNA to <br />
#' _EGAR00001456313_GBM9_R1_SC93.RNA (recurrent tumor)
#' * Blood (normal) samples: _EGAR00001456508_GBM9.DNA_B
#' * Bulk samples (DNA WES): <br />
#' _EGAR00001456467_GBM9_1.DNA (initial tumor), <br />
#' _EGAR00001456468_GBM9_2.DNA (initial tumor), <br />
#' _EGAR00001456469_GBM9_R1.DNA (recurrent tumor), <br />
#' _EGAR00001456470_GBM9_R2.DNA (recurrent tumor)
#'
#' \code{GBM9_postproc} &nbsp; An object of class \code{S4}. Access internal
#' objects from slots using @.
#'
#' \code{GBM9_postproc@annovar.qc} &nbsp; A 23 x 39 dataframe containing the
#' post-processed functional annotations of the somatic variants. For details
#' regarding the column headers, we refer readers to the extensive
#' ANNOVAR [documentation](https://annovar.openbioinformatics.org/) and the
#' original research article \insertCite{Wang2010}{Canopy2}. Only those
#' annotations corresponding to mutations with at least one non-zero count
#' across all cells were retained.
#'
#' \code{GBM9_postproc@featurecounts.qc} &nbsp; A 63677 x 127 matrix of
#' post-processed mapped reads from the single-cell gene expression data.
#' The row names denote the gene names. The column names from the third column
#' onward denote the scRNA-seq samples from the gene expression data for patient
#' GBM9. Only those cells that have at least one non-zero count across genes
#' were retained.
#'
#' \code{GBM9_postproc@param.est} &nbsp; A list of 5 elements: \code{alpha},
#' \code{beta}, \code{scale}, \code{id.g}, and \code{pct.estimable}. The first
#' three are of length equal to the number of mutations where \code{alpha} denotes the
#' mutation-specific gene activation rates, \code{beta} denotes the mutation-specific
#' gene deactivation rates, and \code{scale} denotes the mutation-specific
#' transcription rates. The ids in \code{id.g} denote the positions corresponding
#' to the mutations that were estimable (returning non-NA or non-negative values),
#' and \code{pct.estimable} denotes the percentage that were estimable.
#'
#' \code{GBM9_postproc@Rb} &nbsp; A 18 x 4 matrix of post-processed alternative
#' read counts from the bulk DNA WES data. The row names denote the chromosomal
#' position of the point mutations. The column names denote the bulk samples
#' according to the above labels from the SNV callset of patient GBM9.
#' Post-processing involved removing mutations (rows) with zero alternative
#' read counts across all cells from the single-cell data.
#'
#' \code{GBM9_postproc@Rs} &nbsp; A 18 x 125 matrix of post-processed alternative
#' read counts from the scRNA-seq data. The row names denote the chromosomal
#' position of the point mutations. The column names denote the scRNA-seq samples
#' according to the above labels from the SNV callset of patient GBM9.
#' Post-processing involved i) removing mutations (rows) with zero alternative
#' read counts across all cells and ii) removing cells (columns) with zero total
#' (benign + mutated) read counts across all mutations.
#'
#' \code{GBM9_postproc@Xb} &nbsp; A 18 x 4 matrix of post-processed total read
#' counts (benign + mutated) from the bulk DNA WES data. The row names denote
#' the chromosomal position of the point mutations. The column names denote the
#' bulk samples according to the above labels from the SNV callset of patient
#' GBM9. Post-processing involved removing mutations (rows) with zero
#' alternative read counts across all cells from the single-cell data.
#'
#' \code{GBM9_postproc@Xs} &nbsp; A 18 x 125 matrix of post-processed total read
#' counts (benign + mutated) from the scRNA-seq data. The row names denote the
#' chromosomal position of the point mutations. The column names denote the
#' scRNA-seq samples according to the above labels from the SNV callset of
#' patient GBM9. Post-processing involved i) removing mutations (rows) with zero
#' alternative read counts across all cells and ii) removing cells (columns)
#' with zero total (benign + mutated) read counts across all mutations.
#'
#'@details
#'
#'  | ***Patient*** | ***Subtype*** | ***# bulk DNA WES*** | ***# bulk RNA-seq*** | ***# scRNA-seq (post-processed)***|
#'  |-----	|-----	|-----	|-----	|-----	|
#'  |-----	|-----	|-----	|-----	|-----	|
#'  | GBM9	| Blood                     	| 1              	| 0              	| 0                	|
#'  |       | GBM9_1 (initial tumor)    	| 1              	| 1              	| 29 (25)          	|
#'  |       | GBM9_2 (initial tumor)    	| 1              	| 1              	| 60 (58)          	|
#'  |       | GBM9_R1 (recurrent tumor) 	| 1              	| 1              	| 44 (43)          	|
#'  |       | GMB9_R2 (recurrent tumor) 	| 1              	| 1              	| 0                	|
#'  |       | (total)                   	| 5              	| 4              	| 133 (125)        	|
#'
#' @md
#'
#' @import Rdpack
#' @references{
#'    \insertAllCited{}
#' }
#'
"GBM9_postproc"
#'
#' Pre-processed glioblastoma (GBM) data from patient GBM10
#'
#' Pre-processed single-cell RNA sequencing (scRNA-seq), bulk RNA sequencing
#' (bulk RNA-seq) and bulk DNA whole exome sequencing data (bulk DNA WES) from
#' one patient (GBM10) selected from 52 individuals with glioblastoma (GBM).
#' The authors characterized genomic and expression profiles across 127
#' multisector or longitudinal specimens with matched tumor/normal samples
#' \insertCite{Lee2017}{Canopy2}. The single-nucleotide variant (SNV) pipeline
#' utilized to pre-process the data is described in the main text and supplement.
#' Tumor samples were acquired based on 4 categories: 1) locally adjacent
#' tumors, 2) multifocal/multicentric tumors, 3) 5-ALA (+/-) tumors and
#' 4) Longitudinal tumors. <br />
#' Accession number: EGAS00001001880
#'
#' @format
#'
#' Labels:
#' * Bulk samples (RNA-seq): <br />
#' _EGAR00001456512_GBM10_1.RNA, <br />
#' _EGAR00001456513_GBM10_2.RNA, <br />
#' _EGAR00001456514_GBM10_3.RNA
#' * Single-cell samples (RNA-seq): <br />
#' _EGAR00001455959_GBM10_1_SC06.RNA to _EGAR00001456043_GBM10_2_SC94.RNA
#' * Blood (normal) samples: <br />
#' _EGAR00001456508_GBM10.DNA_B
#' * Bulk samples (DNA WES): <br />
#' _EGAR00001456471_GBM10_1.DNA, <br />
#' _EGAR00001456473_GBM10_2.DNA <br />
#' _EGAR00001456472_GBM10_3.DNA
#'
#' @format
#'
#' \code{GBM10_preproc} &nbsp; An object of class \code{S4}. Access internal
#' objects from slots using @.
#'
#' \code{GBM10_preproc@alt.qc} &nbsp; A 47 x 84 matrix of pre-processed alternative
#' read counts from the bulk and single-cell data for GBM10. The row names denote
#' the chromosomal position of the point mutations. The column names denote the
#' bulk RNA-seq, scRNA-seq, bulk DNA WES, and blood samples according to the
#' above labels from the SNV callset of patient GBM10.
#'
#' \code{GBM10_preproc@annovar} &nbsp; A 47 x 39 dataframe containing the
#' functional annotation of the somatic variants for GBM10. For details regarding
#' the column headers, we refer readers to the extensive
#' ANNOVAR [documentation](https://annovar.openbioinformatics.org/) and the
#' original research article \insertCite{Wang2010}{Canopy2}.
#'
#' \code{GBM10_preproc@featurecounts.qc} &nbsp; A 63677 x 77 matrix of
#' pre-processed mapped reads from single-cell gene expression data. The row
#' names denote the gene names. The column names from the third column onward
#' denote the scRNA-seq samples from the gene expression data for patient GBM10.
#'
#' \code{GBM10_preproc@mut} &nbsp; A 47 x 5 dataframe containing the chromosomal
#' positions of point mutations for GBM10. The columns consist of: the chromosome
#' number \code{Chr}, the starting position of the point mutation \code{Start},
#' the ending position of the point mutation \code{End} (coincides with
#' \code{Start} because the mutation is a single nucleotide variant), the
#' reference allele \code{Ref}, and the alternate allele \code{Alt}.
#'
#' \code{GBM10_preproc@ref.qc} &nbsp; A 47 x 84 matrix of pre-processed
#' reference read counts from the bulk and single-cell data for GBM10. The
#' row names denote the chromosomal position of the point mutations. The
#' column names denote the bulk RNA-seq, scRNA-seq, bulk DNA WES, and blood
#' samples from the total read counts of patient GBM10.
#'
#' \code{GBM10_preproc@summary.qc} &nbsp; A 14 x 77 matrix of the
#' \code{featureCounts} output summary that reports assignment of alignments
#' to genomic features. The row names denote the assignment of the alignments.
#' The column names denote the scRNA-seq samples from the gene expression data
#' for patient GBM10.
#'
#'@details
#'
#'  | ***Patient*** | ***Subtype*** | ***# bulk DNA WES*** | ***# bulk RNA-seq*** | ***# scRNA-seq (pre-processed)***|
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
"GBM10_preproc"
#'
#' Post-processed breast cancer (BC) data from patient GBM10
#'
#' Some further data munging of GBM10_preproc to get GBM10_postproc, which contains
#' the finalized read counts for Canopy2 functions and vignette.<br />
#' Accession number: EGAS00001001880
#'
#' @format
#'
#' Labels:
#' * Bulk samples (RNA-seq): <br />
#' _EGAR00001456512_GBM10_1.RNA, <br />
#' _EGAR00001456513_GBM10_2.RNA, <br />
#' _EGAR00001456514_GBM10_3.RNA
#' * Single-cell samples (RNA-seq): <br />
#' _EGAR00001455959_GBM10_1_SC06.RNA to _EGAR00001456043_GBM10_2_SC94.RNA
#' * Blood (normal) samples: <br />
#' _EGAR00001456508_GBM10.DNA_B
#' * Bulk samples (DNA WES): <br />
#' _EGAR00001456471_GBM10_1.DNA, <br />
#' _EGAR00001456473_GBM10_2.DNA, <br />
#' _EGAR00001456472_GBM10_3.DNA
#'
#'
#' \code{GBM10_postproc} &nbsp; An object of class \code{S4}. Access internal
#' objects from slots using @.
#'
#' \code{GBM10_postproc@annovar.qc} &nbsp; A 8 x 39 dataframe containing the
#' post-processed functional annotations of the somatic variants. For details
#' regarding the column headers, we refer readers to the extensive
#' ANNOVAR [documentation](https://annovar.openbioinformatics.org/) and the
#' original research article \insertCite{Wang2010}{Canopy2}. Only those
#' annotations corresponding to mutations with at least one non-zero count
#' across all cells were retained.
#'
#' \code{GBM10_postproc@featurecounts.qc} &nbsp; A 63677 x 73 matrix of
#' post-processed mapped reads from the single-cell gene expression data.
#' The row names denote the gene names. The column names from the third column
#' onward denote the scRNA-seq samples from the gene expression data for patient
#' GBM10. Only those cells that have at least one non-zero count across genes
#' were retained.
#'
#' \code{GBM10_postproc@param.est} &nbsp; A list of 5 elements: \code{alpha},
#' \code{beta}, \code{scale}, \code{id.g}, and \code{pct.estimable}. The first
#' three are of length equal to the number of mutations where \code{alpha} denotes the
#' mutation-specific gene activation rates, \code{beta} denotes the mutation-specific
#' gene deactivation rates, and scale denotes the mutation-specific
#' transcription rates. The ids in \code{id.g} denote the positions corresponding
#' to the mutations that were estimable (returning non-NA or non-negative values),
#' and \code{pct.estimable} denotes the percentage that were estimable.
#'
#' \code{GBM10_postproc@Rb} &nbsp; A 8 x 3 matrix of post-processed alternative
#' read counts from the bulk DNA WES data. The row names denote the chromosomal
#' position of the point mutations. The column names denote the bulk samples
#' according to the above labels from the SNV callset of patient GBM10.
#' Post-processing involved removing mutations (rows) with zero alternative
#' read counts across all cells from the single-cell data.
#'
#' \code{GBM10_postproc@Rs} &nbsp; A 8 x 71 matrix of post-processed alternative
#' read counts from the scRNA-seq data. The row names denote the chromosomal
#' position of the point mutations. The column names denote the scRNA-seq samples
#' according to the above labels from the SNV callset of patient GBM10.
#' Post-processing involved i) removing mutations (rows) with zero alternative
#' read counts across all cells and ii) removing cells (columns) with zero total
#' (benign + mutated) read counts across all mutations.
#'
#' \code{GBM10_postproc@Xb} &nbsp; A 8 x 3 matrix of post-processed total read
#' counts (benign + mutated) from the bulk DNA WES data. The row names denote
#' the chromosomal position of the point mutations. The column names denote the
#' bulk samples according to the above labels from the SNV callset of patient
#' GBM10. Post-processing involved removing mutations (rows) with zero
#' alternative read counts across all cells from the single-cell data.
#'
#' \code{GBM10_postproc@Xs} &nbsp; A 8 x 71 matrix of post-processed total read
#' counts (benign + mutated) from the scRNA-seq data. The row names denote the
#' chromosomal position of the point mutations. The column names denote the
#' scRNA-seq samples according to the above labels from the SNV callset of
#' patient GBM10. Post-processing involved i) removing mutations (rows) with zero
#' alternative read counts across all cells and ii) removing cells (columns)
#' with zero total (benign + mutated) read counts across all mutations.
#'
#'@details
#'
#'  | ***Patient*** | ***Subtype*** | ***# bulk DNA WES*** | ***# bulk RNA-seq*** | ***# scRNA-seq (post-processed)***|
#'  |-----	|-----	|-----	|-----	|-----	|
#'  |-----	|-----	|-----	|-----	|-----	|
#'  | GBM10 | Blood                     	| 1              	| 0              	| 0                	|
#'  |       | GBM10_1                   	| 1              	| 1              	| 31 (28)          	|
#'  |       | GBM10_2                   	| 1              	| 1              	| 54 (43)          	|
#'  |       | GBM10_3                   	| 1              	| 1              	| 0                	|
#'  |       | (total)                   	| 4              	| 3              	| 85 (71)          	|
#'
#' @md
#'
#' @import Rdpack
#' @references{
#'    \insertAllCited{}
#' }
#'
"GBM10_postproc"
