# MISTy companion data
# Copyright (c) 2020 Jovan Tanevski, Attila Gabor <attila.gabor@uni-heidelberg.de>

#' Synthetic benchmark data for MISTy
#'
#' Data generated from 10 random layouts of four cell types and empty space on
#' 100-by-100 gridby simulating a two-dimensional cellular automata model that
#' focuses on signaling events. Cell growth, division, motility and death are
#' neglected. The intracellular processes involve two layers, first the ligand
#' activation of signaling hubs and ligand production/secretion regulated by
#' proteins. The model simulates the production, diffusion, degradation and
#' interactions of 11 molecular species. Ligands are produced in each cell-type
#' based on the activity level of their production nodes and then freely diffuse,
#' degrade or interact with other cells on the grid. Other molecular species
#' involved in signaling are localised in the intracellular space and their
#' activity depends on ligand binding and intracellular wiring.
#'
#' @format A named \code{list} of length 10. Each list item is a \code{tibble} that
#' corresponds to a simulation of one random layout with information about each
#' cell in rows described by the following 14 variables:
#' \describe{
#'  \item{row, col}{location of the cell on the grid}
#'  \item{ligA, ligB, ligC, ligD}{expression of ligands}
#'  \item{protE, protF}{expression of intracellular proteins}
#'  \item{prodA, prodB, prodC, prodD}{expression of regulatory proteins}
#'  \item{type}{cell type id}
#' }
#'
#'
#' @source \url{https://github.com/saezlab/misty_pipelines/}
"synthetic"
