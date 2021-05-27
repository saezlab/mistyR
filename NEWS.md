# mistyR 1.0.2

- Bugfix: models built with different parameters stored and retrieved from the same cache file.
- Avoid calls to os-dependent file.info in tests. 

# mistyR 1.0.1

- Bugfix: passing arguments to ranger.
- Warnings on clearing nonexistent cache folders and tests of performance.
- Increased test coverage.

# mistyR 1.0.0

- Release version for Bioconductor 3.13.

# mistyR 0.99.11

- Fixed a bug in Nystrom approximation for creating paraview.
- Added a suite of tests with high coverage.
- Cleaner cache control.

# mistyR 0.99.9

## Revisions requested by Bioconductor

- Removed magrittr from dependencies. Reexported pipe operator. Removed from examples.
- Caching turned off by default.
- Added parameter for appending performance and coefficient files in run_misty.
- Internal functions in views.R are now explicit.
- Removed alternatives and additional examples for function use from documentation.
- Removed redundant messages, escalated to warnings where required.
- run_misty cleans up empty cache directories.
- Vignette compatibility with the new release of SpatialExperiment.
- Remove Seurat vignette from package due to missing hdf5r binary for R 4.1 on CRAN for MacOS.

## Other changes

- README.md figures moved to the cloud.
- All passed paths and cache location are normalized.
- Changes in the mistyR vignette to reflect changes to insilico evaluation from the paper.

# mistyR 0.99.0

-   Version with vignettes ready to submit to Bioconductor.

# mistyR 0.1.0 (MISTy)

-   Initial beta release of mistyR (named as MISTy) with function documentation.