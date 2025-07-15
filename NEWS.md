# PopComm 0.1.1.0 (2025-07-15)

## Changes in circle_plot()

- Updated `circle_plot()` to support compatibility with `igraph` version 2.1.5 and above.
- For `igraph <= 2.1.4`, self-loop angles are calculated manually to radiate outward based on node positions.
- For `igraph >= 2.1.5`, loop angles are handled automatically by `igraph`; added logic to suppress or delete self-loops when disabled.
- Added internal version check using `utils::packageVersion()` to differentiate logic paths and ensure backward compatibility.
- Removed direct assignment to `E(g)$loop.angle` for new `igraph` versions, which now throw errors with incorrect defaults.

## Miscellaneous

- Added `@importFrom utils packageVersion` directive to support version checking.
- Improved robustness of layout-based self-loop rendering.

