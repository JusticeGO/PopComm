#' Ligand-Receptor Pair Database
#'
#' A comprehensive database of human ligand-receptor pairs with gene/protein identifiers
#' and supporting evidence from literature. Data imported from `human_lr_pair.txt`.
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
#' @source Source from CellTalkDB (\url{http://tcm.zju.edu.cn/celltalkdb/}).
"lr_db"


#' Example Seurat Object for PopComm Package
#'
#' A preprocessed Seurat object containing single-cell RNA sequencing data for demonstration
#' and testing purposes across functions in the PopComm package. The object includes normalized
#' expression data, metadata with sample and cell type annotations, and basic preprocessing
#' steps (e.g., dimensionality reduction).
#'
#' @format A \code{\link[Seurat]{Seurat}} object with the following slots:
#' \describe{
#'   \item{assays}{RNA count matrix normalized by log1p.}
#'   \item{meta.data}{Dataframe containing metadata columns:
#'     \itemize{
#'       \item \code{sample}: Sample identifier (e.g., "sample1", "sample2").
#'       \item \code{broad.cell.type}: Cell type annotations (e.g., "Cardiac", "Fibroblast").
#'     }
#'   }
#'   \item{reductions}{PCA and UMAP embeddings for visualization.}
#' }
#' @source Subset of snRNA-seq data analyzed in our study,
#' processed via standard workflow (see manuscript for details).
#' Full dataset available at \href{https://zenodo.org/record/XXXXX}{Zenodo}
#' (DOI: 10.XXXX/zenodo.XXXXX).
#' @keywords internal
#' @export
load_example_seurat <- function() {
  url <- "https://sandbox.zenodo.org/records/167208/files/example_seurat_obj.rds"
  temp_file <- tempfile(fileext = ".rds")

  old_timeout <- getOption("timeout")
  options(timeout = 180)
  on.exit(options(timeout = old_timeout))

  tryCatch(
    {
      utils::download.file(url, temp_file, quiet = TRUE, mode = "wb")
      readRDS(temp_file)
    },
    error = function(e) {
      message("Download failed, please check the network or URL validity. Error message:")
      message(e)
      return(NULL)
    }
  )
}
