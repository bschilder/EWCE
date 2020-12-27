#' @details
#' The EWCE package is designed to facilitate expression weighted cell type enrichment analysis
#' as described in our Frontiers in \href{https://doi.org/10.3389/fnins.2016.00016}{Neuroscience paper}.
#'
#' The package was originally designed to work with the single-cell cortical transcriptome data
#'  from the \href{http://linnarssonlab.org}{Linnaerson Lab}, available \href{http://mousebrain.org/}{here}.
#'
#' Using this package it is possible to read in any single cell transcriptome data,
#'  provided that you have a cell by gene expression matrix (with each cell as a seperate column)
#'  and a seperate annotation dataframe, with a row for each cell.
#' The EWCE process involves testing for whether the genes in a target list have higher levels
#' of expression in a given cell type than can reasonably be expected by chance.
#'
#' The probability distribution for this is estimated by randomly generating gene lists of equal length
#' from a set of background genes.
#' The \pkg{EWCE} method can be applied to any gene list. In the paper we reported its application
#'  to genetic and transcriptomic datasets, and in this vignette we detail how this can be done.
#' Note that throughout this vignette we use the terms ‘cell type’ and ‘sub-cell type’ to refer to
#' two levels of annotation of what a cell type is. This is described in further detail in our \href{https://doi.org/10.3389/fnins.2016.00016}{Neuroscience paper},
#'  but relates to the two levels of annotation provided in the \href{https://pubmed.ncbi.nlm.nih.gov/25700174/}{Linnarsson dataset}.
#'   In this dataset a cell is described as having a cell type (i.e. ‘Interneuron’) and
#'   subcell type (i.e. ‘Int11’ a.k.a Interneuron type 11).
"_PACKAGE"
