# MISTy - **M**ultiview **I**ntercellular **S**pa**T**ial modeling framework <img src="man/figures/logo.png" align="right" height="139">

<!-- badges: start -->
<!-- badges: end -->

## Overview

The advancement of technologies for measurement of highly multiplexed spatial data require the development of scalable methods that can leverage the availability of the spatial context. Multiview Intercellular SpaTial modeling framework (MISTy) is an explainable machine learning framework for knowledge extraction and analysis of single-cell, highly multiplexed, spatially resolved data.

<img src="man/figures/graphical_abstract.png" align="center" width="800">

MISTy facilitates an in-depth understanding of marker interactions by profiling the intra- and intercellular relationships. MISTy is a flexible framework able to process a custom number of views. Each of these views can describe a different spatial context, i.e., define a relationship among the observed expressions of the markers, such as intracellular regulation or paracrine regulation. However, the views can also capture cell-type specific relationships, capture relations between functional footprints or focus on relations between different anatomical regions. Each MISTy view is considered as a potential source of variability in the measured marker expressions. Each MISTy view is then analyzed for its contribution to the total expression of each marker and is explained in terms of the interactions with other measurements that led to the observed contribution. Our approach is modular, easily parallelizable and thus scalable to samples with millions of cells and thousands of measured markers.


## System Requirements

MISTy requires a standard configuration and enough RAM to store the analyzed dataset and to support in-memory operations.

The package requires installed R version 3.6.2 or higher. This package has been tested on a macOS (10.15.4). The package should be compatible with Windows, Linux and maxOS operating systems.


## Installation

Install from GitHub using devtools:

```r
# install.packages("devtools")
devtools::install_github("saezlab/misty")

```

MISTy is dependent on the following packages that are available from CRAN:

MASS (>= 7.3.51.5),
dplyr (>= 0.8.5),
purrr (>= 0.3.4),
furrr (>= 0.1.0),
readr (>= 1.3.1),
stringr (>= 1.4.0),
tibble (>= 3.0.1),
caret (>= 6.0.86),
randomForest (>= 4.6.14),
deldir (>= 0.1.25),
distances (>= 0.1.8),
digest (>= 0.6.25),
rlist (>= 0.4.6.1),
assertthat (>= 0.2.1),
magrittr (>= 1.5),
rlang (>= 0.4.5)

and suggests ranger (>= 0.12.1) and knitr (>= 1.28)

The installation of MISTy without any additional dependencies should take only seconds.

## Usage

Example pipelines and synthetic data for running MISTy are available from [this repository](https://github.com/saezlab/misty_pipelines/). To run MISTy on the provided synthetic data run the script synthetic_pipeline.R.

## Citation
If you use MISTy for your research please cite the [following publication](https://doi.org/10.1101/2020.05.08.084145): 

> Jovan Tanevski, Attila Gabor, Ricardo Omar Ramirez Flores, Denis Schapiro, Julio Saez-Rodriguez. Explainable multi-view framework for dissecting inter-cellular signaling from highly multiplexed spatial data. bioRxiv 2020.05.08.084145 doi: [10.1101/2020.05.08.084145](https://doi.org/10.1101/2020.05.08.084145) (2020).
