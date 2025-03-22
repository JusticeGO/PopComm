# data-raw/generate_data.R
# lr_database
lr_db <- read.delim("data-raw/human_lr_pair.txt")
usethis::use_data(lr_db, overwrite = TRUE)

# Example
# seurat_object <- readRDS("data-raw/example_seurat_obj.rds")
# usethis::use_data(seurat_object)


# Example filter_lr_eg
filtered_lr_eg <- readRDS("data-raw/filtered_lr_eg.rds")
usethis::use_data(filtered_lr_eg, overwrite = TRUE)

# Example lr_scores_eg
lr_scores_eg <- readRDS("data-raw/lr_scores_eg.rds")
usethis::use_data(lr_scores_eg, overwrite = TRUE)

# Example metadata_eg
metadata_eg <- read.csv("data-raw/metadata_eg.tsv", sep = "\t")
usethis::use_data(metadata_eg, overwrite = TRUE)
