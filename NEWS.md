# PopComm Change Logs

## PopComm 0.1.1.0 (2025-07-15)

### Changes in circle_plot()

- Updated `circle_plot()` to support compatibility with `igraph` version 2.1.5 and above.
- For `igraph <= 2.1.4`, self-loop angles are calculated manually to radiate outward based on node positions.
- For `igraph >= 2.1.5`, loop angles are handled automatically by `igraph`; added logic to suppress or delete self-loops when disabled.
- Added internal version check using `utils::packageVersion()` to differentiate logic paths and ensure backward compatibility.
- Removed direct assignment to `E(g)$loop.angle` for new `igraph` versions, which now throw errors with incorrect defaults.

### Miscellaneous

- Added `@importFrom utils packageVersion` directive to support version checking.
- Improved robustness of layout-based self-loop rendering.


## PopComm 0.1.2.0 (2025-07-29)

### Changes in filter_lr_single(), filter_lr_all(), one_step_single(), one_step_all()

- Added new parameters: `min_r2` and `min_fstat` to filter models by minimum R-squared and F-statistic thresholds.
- Output now includes r_squared (`r2`) and f_statistic (`fstat`) for each linear regression model.
- Models with R² below `min_r2` or F-statistic below `min_fstat` are excluded from results.
- Enhances model screening based on goodness-of-fit and overall regression significance.


## PopComm 1.0.0.0 (2025-08-12)

### Changes in filter_lr_single(), filter_lr_all(), score_lr_single(), single_lr_all(), one_step_single(), one_step_all()

- Added support for average expression matrices. These functions now accept either Seurat objects or average expression matrices (a numeric matrix of gene expression where columns encode cell types and samples, and rows are genes).
- No changes to existing Seurat-based workflows.
