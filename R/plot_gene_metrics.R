
#' plot_gene_metrics
#'
#' Convenience functiont to quickly plot multiple aspects of a CellTypeData (CTD) object.
#' @param ctd CellTypeData object produced by \pkg{EWCE}.
#' @param level Which annotation level to plot.
#' @param metric "mean_exp", "specificity", or "specificity_quantiles".
#' @param n_genes Limit plot to just the top N genes (according to \code{metric}).
#' @param title Plot title.
#' @param interactive Make the plot interactive with \pkg{plotly}.
#' @param remove_numeric_groups Easy way to remove unlabeled levels (i.e. clusters only identified by numbers).
#' @examples
#' data("ctd")
#' gp <- plot_gene_metrics(ctd = ctd, metric = "specificity", level = 2, n_genes = 10)
#' @export
#' @import ggplot2
#' @import dplyr
plot_gene_metrics <- function(ctd,
                              level=2,
                              metric="specificity",
                              n_genes=NULL,
                              title=metric,
                              interactive=F,
                              remove_numeric_groups=F,
                              show_plot=T){
    top_markers <- find_celltype_markers(ctd = ctd, level=level)
    if(!is.null(n_genes)) top_markers <- top_markers[1:n_genes]
    marker_df <- ctd[[level]][[metric]][top_markers,]
    marker_df_select <- suppressWarnings(data.frame(marker_df[,is.na(as.numeric(colnames(marker_df)))]))
    exp <- reshape2::melt(cbind(marker_df_select, Gene=top_markers),
                          id.vars="Gene",
                          variable.name="Cell_type",
                          value.name = metric)
    exp <- exp %>%
        dplyr::mutate(factor(Cell_type,
                             levels = sort(unique(Cell_type)), ordered = T ))

    # library(ggplot2)
    gp <- ggplot(data = exp, aes_string(x="Cell_type", y=metric, fill="Cell_type")) +
        geom_bar(stat="identity", alpha=.75, show.legend = F) +
        facet_grid(Gene~.) +
        labs(title=title) +
        scale_y_continuous(n.breaks = 3) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    if(interactive) gp <- plotly::ggplotly(gp)
    if(show_plot) print(gp)
    return(gp)
}
