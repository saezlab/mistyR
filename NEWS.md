# mistyR 0.99.8

## Revisions requested by Bioconductor

-   Removed magrittr from dependencies. Reexported pipe operator. Removed from examples.
-   Caching turned off by default.
-   Added parameter for appending performance and coefficient files in run_misty.
-   Internal functions in views.R are now explicit.
-   Removed alternatives and additional examples for function use from documentation.
-   Removed redundant messages, escalated to warnings where required.

## Minor changes

-   All passed paths and cache location are normalized.
-   run_misty cleans up empty cache directories.

# mistyR 0.99.7

-   Vignette compatibility with the new release of SpatialExperiment.

# mistyR 0.99.6

-   Minor revisions for bioc release.

# mistyR 0.99.5

-   Remove Seurat vignette from package due to missing hdf5r binary for R 4.1 on CRAN for MacOS.

# mistyR 0.99.4

-   Triggering bioc rebuild.

# mistyR 0.99.3

-   Changes in the mistyR vignette to reflect changes to insilico evaluation from the paper.

# mistyR 0.99.2

-   Fix vignette build related to older version of SpatialExperiment.

# mistyR 0.99.1

-   Fix vignette build issues on some systems.

# mistyR 0.99.0

-   Version with vignettes ready to submit to Bioconductor.

# mistyR 0.1.0 (MISTy)

-   Initial beta release of mistyR (named as MISTy) with function documentation.