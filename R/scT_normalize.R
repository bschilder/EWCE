
#' Normalize expression matrix
#'
#' Normalize expression matrix by accounting for library size.
#' Uses \pkg{sctransform}.
#' @export
#' @import sctransform
#' @import Matrix
scT_normalize <- function(exprs){
    scT <- sctransform::vst(exprs,
                            return_cell_attr = T)
    exp_scT <- sctransform::correct_counts(scT, exprs) # UMI_corrected
    exp_scT_normed <- Matrix::t(Matrix::t(exp_scT)*(1/Matrix::colSums(exp_scT)))
    return(exp_scT_normed)
}
