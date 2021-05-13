# MISTy - **M**ultiview **I**ntercellular **S**pa**T**ial modeling framework <img src="man/figures/logo.png" align="right" height="139">

<!-- badges: start -->

| branch | status                                                       |
| :----- | :----------------------------------------------------------- |
| master | [![Build Status](https://travis-ci.org/saezlab/mistyR.svg?branch=master)](https://travis-ci.org/saezlab/mistyR) |
| devel  | [![Build Status](https://travis-ci.org/saezlab/mistyR.svg?branch=devel)](https://travis-ci.org/saezlab/mistyR) |

<!-- badges: end -->

## Overview

The advancement of technologies for measurement of highly multiplexed spatial data require the development of scalable methods that can leverage the availability of the spatial context. Multiview Intercellular SpaTial modeling framework (MISTy) is an explainable machine learning framework for knowledge extraction and analysis of single-cell, highly multiplexed, spatially resolved data.

<img src="man/figures/graphical_abstract.png" align="center" width="800">

MISTy facilitates an in-depth understanding of marker interactions by profiling the intra- and intercellular relationships. MISTy is a flexible framework able to process a custom number of views. Each of these views can describe a different spatial context, i.e., define a relationship among the observed expressions of the markers, such as intracellular regulation or paracrine regulation. However, the views can also capture cell-type specific relationships, capture relations between functional footprints or focus on relations between different anatomical regions. Each MISTy view is considered as a potential source of variability in the measured marker expressions. Each MISTy view is then analyzed for its contribution to the total expression of each marker and is explained in terms of the interactions with other measurements that led to the observed contribution. Our approach is modular, easily parallelizable and thus scalable to samples with millions of cells and thousands of measured markers.

**mistyR** is a R package implementing MISTy.


## System Requirements

mistyR requires a standard configuration and enough RAM to store the analyzed dataset and to support in-memory operations.

The package requires installed R version 4.0 or higher. This package has been tested on a macOS Big Sur. The package should be compatible with Windows, Linux and maxOS operating systems.


## Installation

Install from Bioconductor:

```r
if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("mistyR")
```

or install the latest version from GitHub with `remotes`:

```r
# install.packages("remotes")
remotes::install_github("saezlab/mistyR")
```

The installation of mistyR without any additional dependencies and building vignettes should take only seconds.

## Usage

Start by reading `vignette("mistyR")` to learn how to run mistyR. To use mistyR with commonly used objects for spatial omics data see the [articles](https://saezlab.github.io/mistyR/articles/).

Example pipelines and synthetic data for running mistyR are also available from [this repository](https://github.com/saezlab/misty_pipelines). To run mistyR on the provided synthetic data run the script synthetic_pipeline.R.

## Citation
If you use mistyR for your research please cite the [following publication](https://doi.org/10.1101/2020.05.08.084145): 

> Jovan Tanevski, Attila Gabor, Ricardo Omar Ramirez Flores, Denis Schapiro, Julio Saez-Rodriguez (2020). Explainable multi-view framework for dissecting inter-cellular signaling from highly multiplexed spatial data. *bioRxiv*. doi: [10.1101/2020.05.08.084145](https://doi.org/10.1101/2020.05.08.084145)