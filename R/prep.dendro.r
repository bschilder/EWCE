#' prep.dendro
#'
#' \code{prep.dendro} adds a dendrogram to a cts
#'
#' @param ctdIN A single annotLevel of a ctd, i.e. ctd[[1]] (the function is intended to be used via apply)
#' @return A ctd with dendrogram plotting info added
#' @examples
#' ctd = lapply(ctd,EWCE::bin.specificity.into.quantiles,numberOfBins=40)
#' ctd = lapply(ctd,EWCE::prep.dendro)
#' @export
#' @import ggdendro
prep.dendro <- function(ctdIN,
                        verbose=T){
    printer("Preparing dendrograms...",v=verbose)
    binned_file_dist <- stats::dist(t(ctdIN$specificity_quantiles)) # euclidean distances between the rows
    binned_file_dist_hclust <- stats::hclust(binned_file_dist)
    ddata <- ggdendro::dendro_data(binned_file_dist_hclust, type="rectangle")
    ordered_cells <- as.character(ddata$labels$label)
    a1 <- ggplot2::ggplot(ggdendro::segment(ddata)) + ggplot2::geom_segment(ggplot2::aes_string(x="x", y="y", xend="xend", yend="yend")) + ggplot2::coord_flip() +  ggdendro::theme_dendro()
    a1 <- a1 + ggplot2::scale_x_continuous(expand = c(0, 1.3))
    b1 <- ggplot2::ggplot(ggdendro::segment(ddata)) + ggplot2::geom_segment(ggplot2::aes_string(x="x", y="y", xend="xend", yend="yend")) + ggdendro::theme_dendro()
    b1 <- b1 + ggplot2::scale_x_continuous(expand = c(0, 1.3))
    ctdIN$plotting = list()
    ctdIN$plotting$ggdendro_vertical = a1
    ctdIN$plotting$ggdendro_horizontal = b1
    ctdIN$plotting$cell_ordering = ordered_cells
    return(ctdIN)
}



