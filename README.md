# PopComm

<p align="center">
  <a href="https://CRAN.R-project.org/package=PopComm">
    <img src="https://img.shields.io/cran/v/PopComm?color=blue&label=CRAN&logo=R" alt="CRAN Version" />
  </a>
  <a href="https://CRAN.R-project.org/package=PopComm">
    <img src="https://cranlogs.r-pkg.org/badges/grand-total/PopComm?color=green" alt="CRAN Downloads" />
  </a>
  <a href="https://pypi.org/project/PopComm">
    <img src="https://img.shields.io/pypi/v/PopComm?color=blue&label=PyPI&logo=python" alt="PyPI Version" />
  </a>
  <a href="https://pepy.tech/projects/popcomm">
  <img src="https://static.pepy.tech/badge/PopComm" alt="PyPI Downloads" />
  </a>
  <a href="https://github.com/JusticeGO/PopComm">
    <img src="https://img.shields.io/github/v/release/JusticeGO/PopComm?color=blue&label=Github&logo=Github" alt="GitHub Version" />
  </a>
  <a href="https://github.com/JusticeGO/PopComm">
    <img src="https://img.shields.io/github/downloads/JusticeGO/PopComm/total?color=green" alt="GitHub Downloads" />
  </a>
  <img  src="https://img.shields.io/github/license/JusticeGO/PopComm" alt="GitHub License">
</p>

**PopComm** (*Population-Level Cell–Cell Communication Analysis Tools*) is a computational framework designed to quantify **sample-level** ligand–receptor (LR)–mediated cell–cell communication from human large-scale single-cell and single-nucleus RNA-seq datasets. It infers LR interactions between sender/receiver cell types, computes **per-sample communication strength**, and supports association with phenotypes (e.g., aging, disease) as well as comprehensive visualization.

> PopComm is available in **R** and **Python**. It uses a **human** LR database by default, but users can easily swap it for a species-specific database (e.g., mouse) to suit their research needs.


## Overview
<p align="center">
  <img src="https://github.com/JusticeGO/PopComm/blob/main/Overview_PopComm.png" alt="Overview" style="max-width: 100%;">
</p>

PopComm provides three primary modules for population-scale cell–cell communication analysis:

1. **LR Filtering Module**:
   - Identifies expressed ligand–receptor pairs (LR pairs) between sender and receiver cell types, using expression thresholds and cross-sample robustness criteria (e.g., correlation or model-based evidence).
   - Can operate on specific cell-type pairs or all combinations in the dataset.
2. **LR Scoring Module**:
   - Quantifies sample-level communication intensity using a projection-based scoring approach, enabling downstream comparisons across phenotypes (categorical or continuous).
   - Can operate on specific cell-type pairs or all combinations in the dataset.
3. **Visualization Module**:
   - Built-in visualizations to interpret LR interactions at both the LR and sample levels:
     - Circular network diagrams — Summarize cell–cell communication; edges encode interaction count/strength and node size reflects activity.
     - Dot plots — Rank top ligand–receptor pairs per sender–receiver; dot size shows significance and color the correlation/score.
     - Heatmaps — Show sample-wise LR interaction scores, with optional metadata tracks to highlight phenotypes and patterns.
     - PCA projections — Project sample LR scores into 2D to reveal clustering or separation aligned with metadata.
     - Boxplot (group comparison) — Compare LR scores between discrete groups (e.g., High vs Low) and report test statistics.
     - Dotplot (continuous variable) — Scatter LR scores vs a continuous covariate (e.g., IFNscore); add a fit and show R²/p-value.


## Installation

### R

Install the **PopComm** R package from CRAN (stable) or GitHub (development version):

```R
# Stable version (CRAN - Windows/macOS/Linux)
install.packages("PopComm")

# Development version (GitHub - Windows/macOS/Linux)
# install.packages("devtools")      # If not already installed
devtools::install_github("JusticeGO/PopComm")
```

### Python

The Python version of **PopComm** (available on PyPI) mirrors the R workflow (filter → score → visualize):

```bash
# Stable Version (PyPI)
pip install PopComm
```


## Data Requirements

Prepare single-cell or single-nucleus RNA-seq data plus an LR database table. PopComm works in **Seurat mode** or **Matrix mode**.

### Core inputs
- **`rna`**: A **Seurat** object (v4+) **or** a **matrix** of scRNA/snRNA-seq expression.
- **`lr_database`**: A `data.frame` with columns **`ligand_gene_symbol`** and **`receptor_gene_symbol`**.
- **`sender`** and **`receiver`** (Optional): Cell type designated as the **ligand** sender and the **receptor** receiver.
> **Gene identifiers**: Must match the LR database (e.g., **HGNC** for human, **MGI** for mouse); keep naming/case consistent.

### Seurat mode
- **Required metadata columns** in `rna@meta.data`:
  - **`sample_col`**: column name for sample IDs (e.g., `"sample"`).
  - **`cell_type_col`**: column name for cell type (e.g., `"cell_type"`).
- **Used only in Seurat mode**:
  - **`min_cells`**: minimum cells per **(sample × cell type)**.
  - **`min_cell_ratio`**: fraction of cells expressing the ligand/receptor within sender/receiver.

### Matrix mode
- **Matrix columns** represent aggregated expression for **(cell type, sample)** pairs.
- **Column naming**: encode both pieces with a separator, e.g.  
  `Cardiomyocyte--sample_01`, `Fibroblast--sample_02`.
- Set:
  - **`id_sep`**: the separator string (e.g., `"--"`).
  - **`cell_type_col`**, **`sample_col`**: **numeric positions** after splitting the column name by `id_sep`
    (e.g., for `"Cardiomyocyte--sample_01"`, use `cell_type_col=1`, `sample_col=2`).
- **Note**: `min_cells` and `min_cell_ratio` apply to Seurat mode, not Matrix mode.

### Filtering and statistical criteria (set in the function call)
- **`min_samples`**: minimum valid samples to proceed.  
- **`min_sample_ratio`**: fraction of samples where both ligand and receptor are expressed.  
- **`cor_method`**: `"spearman"` (default), `"pearson"`, or `"kendall"`.  
- **`adjust_method`**: multiple-testing correction (default `"BH"`).  
- **`min_adjust_p`**: adjusted p-value cutoff.  
- **`min_cor`**, **`min_r2`**, **`min_fstat`**: optional thresholds on correlation and regression quality.  
- **`num_cores`**: CPU cores for parallelization. Automatically capped at max(1, {system_cores} - 1) if overspecified.
- **`verbose`**: progress messages.

### Metadata for plotting / association
- A **data.frame** with: a **`sample`** column matching the sample IDs used above.
- One or more **variables to analyze** (e.g., age, disease status) used by plotting/stats functions.


## Quick Start (R)
See the [CRAN reference manual](https://CRAN.R-project.org/package=PopComm) for full details.
Example below uses **Perivascular → Endothelial** in **Matrix mode**.

### Data preparation

```r
# Load package and example data
library(PopComm)

data(matrix_object)      # example column name: "Endothelial--sample01"
data(lr_db)
```

### 1) LR filtering

> Filter LR pairs and compute correlations for Perivascular -> Endothelial

```r
filtered_lr <- filter_lr_single(
  rna             = matrix_object,
  sender          = "Perivascular",
  receiver        = "Endothelial",
  lr_database     = lr_db,
  cell_type_col   = 1,      # after splitting colnames by id_sep, cell type is part 1
  sample_col      = 2,      # sample is part 2
  id_sep          = "--",
  min_samples     = 10,
  min_sample_ratio= 0.1,
  cor_method      = "spearman",
  adjust_method   = "BH",
  min_adjust_p    = 0.05,
  min_cor         = 0,
  min_r2          = 0,
  min_fstat       = 0,
  num_cores       = 10,
  verbose         = TRUE
)

if (is.null(filtered_lr) || nrow(filtered_lr) == 0L) {
  stop("No LR pairs passed filters for Perivascular -> Endothelial.")
}

head(filtered_lr)
```

> If you don’t specify a `sender` and `receiver`, use `filter_lr_all()` to analyze **all** sender–receiver combinations.

>  **Output (`filtered_lr`)**
>
> - `ligand`, `receptor` — ligand and receptor gene symbols
> - `cor`, `p_val`, `adjust.p` — correlation and p-values
> - `sender`, `receiver` — sender and receiver cell types
> - `slope`, `intercept`, `r2`, `fstat` — regression metrics

>    *Rows are ordered by ascending `adjust.p` and then descending `cor`.*

### 2) Sample-level scoring

> Compute per-sample LR projection scores using the filtered pairs

```r
lr_scores <- score_lr_single(
  rna            = matrix_object,
  sender         = "Perivascular",
  receiver       = "Endothelial",
  filtered_lr    = filtered_lr,
  cell_type_col  = 1,
  sample_col     = 2,
  id_sep         = "--",
  num_cores      = 10,
  verbose        = TRUE
)

if (is.null(lr_scores) || nrow(lr_scores) == 0L) {
  stop("No LR scores returned; check filters and inputs.")
}

head(lr_scores)
```

> If you used `filter_lr_all()` in the previous step, use `score_lr_all()` here.

>  **Output (`lr_scores`)**
>
> - Inherits key columns from `filtered_lr` (ligand/receptor/sender/receiver, etc.)
> - `sample` — sample identifier
> - `score` — raw projection score
> - `normalized_score` — 0–1 scaled score

>    *Rows are ordered by the `filtered_lr` keys and descending `score`.*


### 3) One-step integrated analysis (optional)

> Run filtering + scoring in one call

```r
res_single <- one_step_single(
  rna             = matrix_object,
  sender          = "Perivascular",
  receiver        = "Endothelial",
  lr_database     = lr_db,
  cell_type_col   = 1,
  sample_col      = 2,
  id_sep          = "--",
  min_samples     = 10,
  min_sample_ratio= 0.1,
  cor_method      = "spearman",
  adjust_method   = "BH",
  min_adjust_p    = 0.05,
  num_cores       = 10,
  verbose         = TRUE
)

if (!is.null(res_single)) {
  lapply(res_single, head)  # res1 (LR stats), res2 (sample scores)
}
```

> For all sender–receiver combinations, use `one_step_all()`.
> Returns two data frames: `res1` (LR statistics) and `res2` (per-sample scores). Returns `NULL` if required cell types/samples are missing or no LR pairs pass thresholds.

------

## Visualization (R)

PopComm includes views such as **circular networks**, **dot plots**, **heatmaps**, **PCA**, and **box/violin** or **scatter/line** for associations.

> **Prepare `metadata`**: a data.frame with a `sample` column matching `lr_scores$sample`, plus variables to analyze.

### Data preparation

```r
# Load package and example data
library(PopComm)

data(filtered_lr_eg)
data(lr_scores_eg)
data(metadata_eg)
```

### 1) LR filtering views

#### `circle_plot`

> Plots a circular LR interaction network with curved directed edges. Nodes are arranged in a circle, and edge widths and colors represent interaction strengths.

```r
p_net <- circle_plot(
  filtered_lr = filtered_lr_eg,
  edge_width = "count",          # <character> one of c("count", "cor")
  node_colors = NULL,
  show_self_interactions = TRUE,
  cutoff = 0
)

print(p_net)
```

#### `dot_plot`

> This function generates a dot plot to visualize LR interaction. Dot sizes are scaled by the correlation coefficient and dot colors represent -log10(adjust.p). The function supports plotting the top interactions per sender-receiver pair or user-specified LR pairs.

```r
p_dot <- dot_plot(
  filtered_lr = filtered_lr_eg,
  top_n = 5,
  axis ="LR-SR",                 # <character> one of c("LR-SR", "SR-LR")
  type_scale = "size",           # <character> one of c("size", "radius")
  selected_LR = NULL
)

print(p_dot)
```

### 2) Sample-level scoring views

#### `heatmap_sample`

> This function generates a heatmap to visualize the LR interaction scores across samples. Rows represent LR pairs and columns represent samples. Optionally, sample metadata can be used to annotate the columns.

```r
p_hm <- heatmap_sample(
  lr_scores = lr_scores_eg,
  metadata = metadata_eg,
  score = "normalized",          # <character> one of c("normalized", "raw")
  selected_sender   = "Endothelial",
  selected_receiver = "Perivascular",
  selected_metadata = c("Sex", "Age_group", "IFN_type"),
  treeheight_row = 50,
  treeheight_col = 50,
  show_LR     = FALSE,
  show_sample = FALSE,
  basic_title = NULL
)

print(p_hm)
```

#### `pca_sample`

> This function performs principal component analysis (PCA) on LR interaction scores across samples, and generates a scatter plot of the first two principal components. Optionally, sample metadata can be used to color the points.

```r
pca_res <- pca_sample(
  lr_scores = lr_scores_eg,
  metadata = metadata_eg,
  selected_sender   = NULL,
  selected_receiver = NULL,
  color_by = "IFN_type",
  n_components = 2
)

print(pca_res$plot)
head(pca_res$df)
```

### 3) Differential interaction analysis

#### `boxplot_lr_group_comparison` (discrete variable)

> This function generates a boxplot comparing LR interaction scores across sample groups with optional significance testing (t-test or Wilcoxon).

```r
p_box <- boxplot_lr_group_comparison(
  lr_scores = lr_scores_eg,
  metadata = metadata_eg,
  ligand = "PSAP",
  receptor = "LRP1",
  sender = "Perivascular",
  receiver = "Fibroblast",
  group_by = "IFN_type",
  score = "normalized",          # <character> one of c("normalized", "raw")
  test = TRUE,
  paired = FALSE,
  test_method = "wilcox.test"    # <character> one of c("wilcox.test", "t.test")
  colors = c("#5fa9d1", "#154778"),
  title = NULL
)

print(p_box$plot)
head(p_box$df)
```

#### `dotplot_lr_continuous_group` (continuous variable)

> This function creates a dotplot (scatter plot) of LR interaction scores against a continuous variable with optional regression line.

```r
p_scatter <- dotplot_lr_continuous_group(
  lr_scores = lr_scores_eg,
  metadata = metadata_eg,
  ligand = "HLA-A",
  receptor = "LILRB2",
  sender = "Lymphoid",
  receiver = "Myeloid",
  group_by = "IFNscore",
  score = "normalized",          # <character> one of c("normalized", "raw")
  point_size = 3,
  point_color = "dodgerblue4",
  add_regression = TRUE,
  title = NULL
)

print(p_scatter$plot)
head(p_scatter$df)
```

------

### Tips

- **Gene IDs** must match the LR database (e.g., HGNC for human); keep case consistent.
- **Parallelism**: `num_cores` automatically capped at max(1, {system_cores} - 1) if overspecified.
- **Matrix mode**: `sample_col`/`cell_type_col` are **numeric positions** after splitting column names by `id_sep`.
- **Seurat mode**: use **column names** (characters) for `sample_col`/`cell_type_col`; `min_cells`/`min_cell_ratio` are Seurat-only.
- For all combinations, use `filter_lr_all()` / `score_lr_all()` / `one_step_all()`.


## LR Database: CellTalkDB

PopComm supports using curated LR tables. This package uses the **CellTalkDB** resource as the default LR database for human (and optionally mouse) LR interactions.

- **How PopComm uses it:**
  - The LR database provides high-confidence ligand–receptor pairs. PopComm leverages these pairs to (i) filter for plausible interactions between sender and receiver cell types and (ii) quantify sample-level communication scores based on those pairs.
- **Citing [CellTalkDB](https://xomics.com.cn/celltalkdb/index.php):**
  - Please cite the CellTalkDB paper in any publication using PopComm that uses this database as the default:
  - **Shao X., Liao J., Li C., et al.** *[CellTalkDB: A manually curated database of ligand-receptor interactions in human and mouse](https://pubmed.ncbi.nlm.nih.gov/33147626/)*. **Briefings in Bioinformatics** (2021).
- **Swapping the LR database:**
   You can use a custom LR table:
  1. Provide a `data.frame`/`table` with at least two columns: `ligand_gene_symbol` and `receptor_gene_symbol`.
  2. Ensure gene identifiers match your Seurat object or expression matrix (e.g., gene symbols).
  3. Pass your table to the argument `lr_database`.


## Performance

- **Parallelization**: Enable multi-core processing where supported by the package functions to speed up LR filtering and scoring on large datasets.


## FAQ

**Q1. What advantages does PopComm offer over single-sample communication tools?**
- It focuses on **population-scale robustness** and **sample-level** scoring, allowing rigorous statistical associations with phenotypes (e.g., age, disease, continuous traits) and powering discoveries that are consistent across many samples.

**Q2. Can I analyze a specific sender–receiver pair only?**
- Yes. Filter and score modules can be restricted to selected cell-type pairs to reduce computation and enhance interpretability.

**Q3. Does PopComm work with non-human data?**
- Yes. Use an LR table appropriate for your organism (or map orthologs to the species in your expression data) and ensure gene identifiers are consistent.


## Cite PopComm & Database Resources

**PopComm:**
 *A Population-Scale Single-Cell Atlas of the Human Heart Reveals Cellular Remodeling and Cell–Cell Communication in Aging and Cardiac Disease*
 (Submission in progress)

**Database Source:**
We thanks the CellTalk Database developed by Professor Xiaohui Fan. If you use the default LR database in PopComm, please also cite the orignal paper of 
[CellTalkDB](https://xomics.com.cn/celltalkdb/index.php): 
 **Shao X., Liao J., Li C., et al.** *[CellTalkDB: A manually curated database of ligand-receptor interactions in human and mouse](https://pubmed.ncbi.nlm.nih.gov/33147626/)*. **Briefings in Bioinformatics** (2021).

 
## License

PopComm is released under the **MIT License**.
 If you distribute an LR database derived from third-party resources (e.g., **CellTalkDB**), please respect their terms and **acknowledge/cite** those resources.
