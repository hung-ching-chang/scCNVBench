# -----------------------------------------------------------------------------
# Source: scDesign3, developed by Dongyuan Song and collaborators.
# Original code: https://github.com/SONGDONGYUAN1994/scDesign3
#
# Reference:
# Song D, et al. (2024). scDesign3 generates realistic in silico data
# for multimodal single-cell and spatial omics.
# Nature Biotechnology, 42, 247-252.
# https://doi.org/10.1038/s41587-023-01772-1
#
# License: MIT. Original copyright (c) 2024 Dongyuan Song.
# Modifications: Adapted for this study; see MODIFICATIONS.md.
# -----------------------------------------------------------------------------

#' A SingleCellExperiment object containing both cell type and pseudotime
#' @format A dataset with 10 rows (genes) and 1289 cols (cells)
#' @return The corresponding SingleCellExperiment object
#' @usage data("example_sce")
"example_sce"

#' A sparse matrix with example data
#' @format A dgCMatrix with 100 rows (genes) and 2087 cols (cells)
#' @return The corresponding sparse matrix object 'example_count'
#' @usage data("example_count")
"example_count"

#' A numeric vector of pseudotime values
#' @format A numeric vector with 2087 values
#' @return The corresponding pseudotime of the 'example_count'
#' @usage data("pseudotime")
"pseudotime"