# data-raw/generate_data.R
# lr_database
lr_db <- read.delim("data-raw/human_lr_pair.txt")
usethis::use_data(lr_db, overwrite = TRUE)

# Example
# seurat_object <- readRDS("data-raw/example_seurat_obj.rds")
# usethis::use_data(seurat_object)
