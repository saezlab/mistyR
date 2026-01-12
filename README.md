# MISTy - **M**ultiview **I**ntercellular **S**pa**T**ial modeling framework <img src="https://www.dropbox.com/s/giluat1vo5wa7jq/misty_badge.png?raw=1" align="right" height="139"/>

<!-- badges: start -->

[![BioC status](http://www.bioconductor.org/shields/build/release/bioc/mistyR.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/mistyR) [![BioC devel status](http://www.bioconductor.org/shields/build/devel/bioc/mistyR.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/mistyR) [![Codecov test coverage](https://codecov.io/gh/saezlab/mistyR/branch/master/graph/badge.svg)](https://codecov.io/gh/saezlab/mistyR?branch=master) [![Docker Cloud Automated build](https://img.shields.io/docker/cloud/automated/tanevski/mistyr)](https://hub.docker.com/r/tanevski/mistyr) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/bioconductor-mistyr/README.html)

<!-- badges: end -->

## Overview

The advancement of technologies for measurement of highly multiplexed spatial data require the development of scalable methods that can leverage the availability of the spatial context. Multiview Intercellular SpaTial modeling framework (MISTy) is an explainable machine learning framework for knowledge extraction and analysis of single-cell, highly multiplexed, spatially resolved data.

<img src="https://www.dropbox.com/s/5wyh520i1wl62ul/graphical_abstract.png?raw=1" align="center" width="800"/>

MISTy facilitates an in-depth understanding of marker interactions by profiling the intra- and intercellular relationships. MISTy is a flexible framework able to process a custom number of views. Each of these views can describe a different spatial context, i.e., define a relationship among the observed expressions of the markers, such as intracellular regulation or paracrine regulation. However, the views can also capture cell-type specific relationships, capture relations between functional footprints or focus on relations between different anatomical regions. Each MISTy view is considered as a potential source of variability in the measured marker expressions. Each MISTy view is then analyzed for its contribution to the total expression of each marker and is explained in terms of the interactions with other measurements that led to the observed contribution. Our approach is modular, easily parallelizable and thus scalable to samples with millions of cells and thousands of measured markers.

**mistyR** is a R package implementing MISTy.

## Recent developments

A **new** and improved version of mistyR is available and actively developed at [jtanevski/mistyR](https://github.com/jtanevski/mistyR). Among other improvements it:

-   *enables advanced relationship signature extraction without intersecting targets*

-   *enables calculation of communities on non-square matrices*

-   *removes file and folder clutter by storing results in database instead of folders*

-   *adds calculation of correlations alongside importances to give an intuition about the sign of the relationship*

This version is also needed to run our new method for analysis of persistent local neighborhoods and downstream translational applications of relationship-based representations - **kasumi**. Check it out at [jtanevski/kasumi](https://github.com/jtanevski/kasumi).

Moving from global robust relationship-based representations, kasumi is a method for the identification of spatially localized neighborhoods in tissues. Following the MISTy framework it extracts local tissue patch-based intra- and intercellular relationships that persistent across samples and conditions. Kasumi learns compressed explainable representations of spatial omics samples while preserving relevant biological signals that are readily deployable for data exploration and hypothesis generation, facilitating translational tasks.

Read more about kasumi in this [publication](https://doi.org/10.1038/s41467-025-59448-0).

## System Requirements

**mistyR** requires a standard configuration and enough RAM to store the analyzed dataset and to support in-memory operations.

The package requires R version 4.0 or higher. This package is developed on macOS Big Sur. The package should be compatible with Windows, Linux and maxOS operating systems.

## Installation

Install from Bioconductor:

``` r
# install.packages("BiocManager")
BiocManager::install("mistyR")
```

You can install the latest stable and development versions from GitHub with `remotes`:

-   stable

``` r
# install.packages("remotes")
remotes::install_github("saezlab/mistyR")
```

-   development

``` r
remotes::install_github("saezlab/mistyR@devel")
```

## Docker

For the released and the latest stable and development versions we also provide [Docker images](https://hub.docker.com/r/tanevski/mistyr) based on the [Rocker project](https://www.rocker-project.org/) - [rocker/r-base](https://github.com/rocker-org/rocker/tree/master/r-base) image.

To create and start a container from the latest docker image and run R in interactive mode:

``` bash
docker run -it tanevski/mistyr:latest
```

## Usage

Start by reading `vignette("mistyR")` to learn how to run **mistyR**. To learn how to use the package with commonly used objects for spatial omics data see the [articles](https://saezlab.github.io/mistyR/articles/).

Example pipelines and synthetic data for **mistyR** are also available from [this repository](https://github.com/saezlab/misty_pipelines). To run **mistyR** on the provided synthetic data run the script *synthetic_pipeline.R*.

## Citation

If you use **mistyR or kasumi** for your research please cite the following publications:

> Jovan Tanevski, Ricardo Omar Ramirez Flores, Attila Gabor, Denis Schapiro, Julio Saez-Rodriguez. Explainable multiview framework for dissecting spatial relationships from highly multiplexed data. Genome Biology 23, 97 (2022). <https://doi.org/10.1186/s13059-022-02663-5>

> Jovan Tanevski, Loan Vulliard, Miguel A. Ibarra-Arellano, Denis Schapiro, Felix J. Hartmann, Julio Saez-Rodriguez. Learning tissue representation by identification of persistent local patterns in spatial omics data. Nat Communications 16, 4071 (2025). <https://doi.org/10.1038/s41467-025-59448-0>
