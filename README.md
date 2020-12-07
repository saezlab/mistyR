# MISTy - **M**ultiview **I**ntercellular **S**pa**T**ial modeling framework <img src="man/figures/logo.png" align="right" height="139">

<!-- badges: start -->
- master [![Build Status](https://travis-ci.org/saezlab/misty.svg?branch=master)](https://travis-ci.org/saezlab/misty)
- devel [![Build Status](https://travis-ci.org/saezlab/misty.svg?branch=devel)](https://travis-ci.org/saezlab/misty)
<!-- badges: end -->

## Overview

The advancement of technologies for measurement of highly multiplexed spatial data require the development of scalable methods that can leverage the availability of the spatial context. Multiview Intercellular SpaTial modeling framework (MISTy) is an explainable machine learning framework for knowledge extraction and analysis of single-cell, highly multiplexed, spatially resolved data.

<img src="man/figures/graphical_abstract.png" align="center" width="800">

MISTy facilitates an in-depth understanding of marker interactions by profiling the intra- and intercellular relationships. MISTy is a flexible framework able to process a custom number of views. Each of these views can describe a different spatial context, i.e., define a relationship among the observed expressions of the markers, such as intracellular regulation or paracrine regulation. However, the views can also capture cell-type specific relationships, capture relations between functional footprints or focus on relations between different anatomical regions. Each MISTy view is considered as a potential source of variability in the measured marker expressions. Each MISTy view is then analyzed for its contribution to the total expression of each marker and is explained in terms of the interactions with other measurements that led to the observed contribution. Our approach is modular, easily parallelizable and thus scalable to samples with millions of cells and thousands of measured markers.


## System Requirements

MISTy requires a standard configuration and enough RAM to store the analyzed dataset and to support in-memory operations.

The package requires installed R version 4.0.3 (Bunny-Wunnies Freak Out) or higher. This package has been tested on a macOS (11.0.1). The package should be compatible with Windows, Linux and maxOS operating systems.


## Installation

Install from GitHub using devtools:

```r
# install.packages("remotes")
remotes::install_github("saezlab/misty")

```

MISTy is dependent on the following packages that are available from CRAN:

assertthat,
caret,
deldir,
digest,
distances,
dplyr,
filelock,
furrr,
future,
ggplot2,
magrittr,
MASS,
purrr,
ranger,
readr,
rlang,
rlist,
stats,
stringr,
tibble,
tidyr

and suggests igraph and knitr.

The installation of MISTy without any additional dependencies should take only seconds.

## Usage

Example pipelines and synthetic data for running MISTy are available from [this repository](https://github.com/saezlab/misty_pipelines/). To run MISTy on the provided synthetic data run the script synthetic_pipeline.R.

## Citation
If you use MISTy for your research please cite the [following publication](https://doi.org/10.1101/2020.05.08.084145): 

> Jovan Tanevski, Attila Gabor, Ricardo Omar Ramirez Flores, Denis Schapiro, Julio Saez-Rodriguez. Explainable multi-view framework for dissecting inter-cellular signaling from highly multiplexed spatial data. bioRxiv 2020.05.08.084145 doi: [10.1101/2020.05.08.084145](https://doi.org/10.1101/2020.05.08.084145) (2020).
