
#' Find the top gene markers for each cell-type
#'
#' @export
find_celltype_markers <- function(ctd,
                                  level=2){
    DF <- ctd[[level]]$specificity
    top_markers <- row.names(DF)[apply(DF,2,which.max)]
    return(top_markers)
}
