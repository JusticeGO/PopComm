PopComm
---

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

**PopComm** (Population-Level Cell-Cell Communication Analysis Tools) is a computational framework designed to quantify cell–cell communication at the population level using large-scale single-cell and single-nucleus RNA-seq datasets.     The package infers ligand–receptor (LR) interactions across different cell types, and scores the strength of these interactions at the sample level. PopComm allows for in-depth analysis of communication patterns in various phenotypic variables and offers tools for LR pair-level, sample-level, and differential interaction analyses, with comprehensive visualization support to aid biological interpretation.

Both **R** and **Python** versions of **PopComm** are available, providing flexibility for users across different programming environments.


## Overview
<p align="center">
  ![Overview](https://github.com/JusticeGO/PopComm/blob/main/Overview_PopComm.png)
</p>

PopComm includes three primary modules for analyzing and visualizing cell–cell communication:

1. **LR Filtering Module**:
   - Identifies expressed ligand–receptor pairs between sender and receiver cell types.
   - Filters pairs based on expression thresholds, correlation across samples, and linear modeling.
   - Applicable to specific cell-type pairs or all combinations in the dataset.
2. **LR Scoring Module**:
   - Quantifies communication intensity at the sample level using a projection-based scoring approach.
   - Enables the comparison of interaction strength across different phenotypes (e.g., age, disease status).
3. **Visualization Module**:
   - Provides functions to visualize ligand-receptor interactions as:
     - Circular networks
     - Dot plots
     - Heatmaps
     - PCA projections
     - box plots
   - Facilitates both cell-type and sample-level interpretation of communication patterns.

## Installation
