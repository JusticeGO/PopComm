#' Ligand-Receptor Pair Database
#'
#' A comprehensive database of human ligand-receptor pairs with gene/protein identifiers
#' and supporting evidence from literature. Data imported from `human_lr_pair.txt`.
#' CellTalkDB: A manually curated database of ligand-receptor interactions in human and mouse
#'
#' @format A data frame with 3,398 rows (pairs) and 10 columns:
#' \describe{
#'   \item{lr_pair}{Character. Unique identifier for ligand-receptor pair, formatted as "LIGAND_RECEPTOR" (e.g., "SEMA3F_PLXNA3")}
#'   \item{ligand_gene_symbol}{Character. Official HGNC symbol of the ligand gene (e.g., "SEMA3F")}
#'   \item{receptor_gene_symbol}{Character. Official HGNC symbol of the receptor gene (e.g., "PLXNA3")}
#'   \item{ligand_gene_id}{Integer. Entrez Gene ID of the ligand gene (NCBI identifier)}
#'   \item{receptor_gene_id}{Integer. Entrez Gene ID of the receptor gene (NCBI identifier)}
#'   \item{ligand_ensembl_protein_id}{Character. Ensembl protein ID of the ligand (e.g., "ENSP00000002829")}
#'   \item{receptor_ensembl_protein_id}{Character. Ensembl protein ID of the receptor (e.g., "ENSP00000358696")}
#'   \item{ligand_ensembl_gene_id}{Character. Ensembl gene ID of the ligand (e.g., "ENSG00000001617")}
#'   \item{receptor_ensembl_gene_id}{Character. Ensembl gene ID of the receptor (e.g., "ENSG00000130827")}
#'   \item{evidence}{Character. PubMed IDs (PMIDs) supporting the interaction, comma-separated (e.g., "15721238")}
#' }
#'
#' @source Source from CellTalkDB (PMID: 33147626).
"lr_db"

#' Example for filtered_lr
"filtered_lr_eg"

#' Example for lr_scores
"lr_scores_eg"

#' Example for metadata
"metadata_eg"

#' Example for matrix object
"matrix_object"

#' #' Example for seurat object
#' "seurat_object"

