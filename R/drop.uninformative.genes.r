#' drop.uninformative.genes
#'
#' \code{drop.uninformative.genes} drops genes from an SCT expression matrix if they do not significantly vary between any celltypes.
#' Makes this decision based on use of an ANOVA (implemented with Limma). If the F-statistic for variation amongst type2 annotations
#' is less than a strict p-threshold, then the gene is dropped.
#'
#' @param exp Expression matrix with gene names as rownames.
#' @param `level2annot Array of cell types, with each sequentially corresponding a column in the expression matrix
#' @param DGE_method Which method to use for the Differential Gene Expression (DGE) step.
#' @param return_sce Whether to return the filtered results
#' as an expression matrix or a \pkg{SingleCellExperiment}.
#' @param min_variance_decile If `min_variance_decile!=NULL`, calculates the variance of the mean gene expression  across `level2annot` (i.e. cell-types),
#' and then removes any genes that are below `min_variance_decile` (on a 0-1 scale).
#' @param adj_pval_thresh Minimum differential expression significance
#' that a gene must demonstrate across `level2annot` (i.e. cell-types).
#' @param drop_nonhuman_genes Whether to drop genes that don't have human orthologues.
#' @param input_species The species that the `exp` dataset comes from
#'  (to be used during the non-orthologues filtering step).
#' @param verbose Whether to print messages (`T` or `F`).
#' @param ... Additional arguments to be passed to the selected DGE method.
#'
#' @return exp Expression matrix with gene names as rownames,
#' or a \pkg{SingleCellExperiment} if \code{return_sce=T}.
#' @examples
#' data("cortex_mrna")
#' cortex_mrna$exp = cortex_mrna$exp[1:300,] # Use only a subset of genes to keep the example quick
#' exp2 = drop.uninformative.genes(exp=cortex_mrna$exp, level2annot=cortex_mrna$annot$level2class)
#' @export
#' @import limma
#' @import stats
#' @source
#' #' \href{https://petehaitch.github.io/BioC2020_DelayedArray_workshop/articles/Effectively_using_the_DelayedArray_framework_for_users.html}{DelayedArray workshop}
drop.uninformative.genes <- function(exp,
                                     level2annot,
                                     DGE_method="limma",
                                     min_variance_decile=NULL,
                                     adj_pval_thresh=0.00001,
                                     drop_nonhuman_genes=F,
                                     input_species=NULL,
                                     sce_save_dir=NULL,
                                     return_sce=F,
                                     verbose=T,
                                     ...){
    DGE_method <- if(is.null(DGE_method)) "" else DGE_method
    #### Convert to SCE format ####
    sce <- ingest_data(obj = exp)

    ### Remove non-expressed genes ####
    printer("+ Removing non-expressed genes...",v=verbose)
    summed <- DelayedArray::rowSums(SummarizedExperiment::assay(sce))
    # Subset the sce object
    sce = sce[summed!=0,]
    printer(paste(nrow(sce)-sum(summed!=0),"/",nrow(sce),
                  "non-expressed genes dropped"), v=verbose)

    #### Remove non-orthologues ####
    if(drop_nonhuman_genes){
        rowDat <- SummarizedExperiment::rowData(sce)
        orths <- convert_orthologues(gene_df=rowDat,
                                     gene_col="rownames",
                                     input_species=input_species,
                                     drop_nonhuman_genes=T,
                                     one_to_one_only=T,
                                     genes_as_rownames=T,
                                     verbose=verbose)
            # Not always sure about how the sce has been named (original species gene names or human orthologues)
           sce <- tryCatch({sce[orths$Gene_orig,]},
                     error=function(e){ sce[orths$Gene,]
                     })
    }

    # #### Normalize ####
    # # Normalize the data before filtering by variance to get a better selection of % genes to remove
    # if(!is.null(normalize_method)){
    #     MAST::zlm( ~ 1)
    # }


    #### Simple variance ####
    if(!is.null(min_variance_decile)){
        # Use variance of mean gene expression across cell types
        # as a fast and simple way to select genes
        sce <- DelayedArray_filter_variance_quantiles(sce = sce,
                                                       level2annot = level2annot,
                                                       n_quantiles = 10,
                                                       min_variance_decile = min_variance_decile,
                                                       verbose = verbose)
    }

    # Make sure the matrix hasn't been converted to characters
    if(class(SummarizedExperiment::assay(sce)[1,1])=="character"){
        # exp <- SummarizedExperiment::assay(sce)
        # exp = as.matrix(exp)
        storage.mode(SummarizedExperiment::assay(sce)) <- "numeric"
    }

    # Run DGE
    start <- Sys.time()
    #### Limma ####
    # Modified original method
    if(tolower(DGE_method)=="limma"){
        eb <- run_limma(sce = sce,
                        level2annot = level2annot,
                        verbose = verbose,
                        ...)
        pF <- stats::p.adjust(eb$F.p.value, method="BH")
        keep_genes <- pF<adj_pval_thresh
        printer(paste(nrow(sce)-sum(keep_genes),"/",nrow(sce),
                      "genes dropped @ DGE adj_pval_thresh <",adj_pval_thresh), v=verbose)
        sce <- sce[keep_genes,]
    }

    #### DESeq2 ####
    if(tolower(DGE_method)=="deseq2"){
        dds_res <- run_DESeq2(sce=sce,
                              level2annot = level2annot,
                              verbose = verbose,
                              ...)
        dds_res <- subset(dds_res, padj<adj_pval_thresh)
        printer(paste(nrow(sce)-nrow(dds_res),"/",nrow(sce),
                      "genes dropped @ DGE adj_pval_thresh <",adj_pval_thresh), v=verbose)
        # Filter original SCE
        sce <- sce[row.names(dds_res), ]
    }

    #### glmGamPoi ####
    ## Removing this option for now until we can figure out how to pass the Travis CI checks,
    ## which are failing when installing the deps for glmGamPo (hdf5).
    # if(tolower(DGE_method)=="glmgampoi"){
    #     # Best for large datasets that can't fit into memory.
    #     messager("DGE:: glmGamPoi...",v=verbose)
    #     sce_de <- run_glmGamPoi_DE(sce,
    #                                level2annot=level2annot,
    #                                adj_pval_thresh=adj_pval_thresh,
    #                                return_as_SCE=T,
    #                                # sce_save_dir=sce_save_dir,
    #                                verbose=verbose)
    #                                # ...)
    #     printer(paste(nrow(sce)-nrow(sce_de),"/",nrow(sce),
    #                   "genes dropped @ DGE adj_pval_thresh <",adj_pval_thresh), v=verbose)
    #     sce <- sce_de
    # }
    # Report time elapsed
    end <- Sys.time()
    print(end-start)
    #### Return results ####
    if(return_sce){
        return(sce)
    }else{
        return(SummarizedExperiment::assay(sce))
    }
}


DelayedArray_filter_variance_quantiles <- function(sce,
                                                   level2annot,
                                                   n_quantiles=10,
                                                   min_variance_decile=.5,
                                                   verbose=T){
    printer("+ Filtering by variance deciles...",v=verbose)
    #### Calculate mean gene expression per cell-type ####
    sce_means <- DelayedArray_grouped_stats(sce = sce,
                                            grouping_var = level2annot,
                                            stat = "mean",
                                            return_sce = F)
    #### Calculate gene variance across cell-types means ####
    gene_variance <- setNames(DelayedMatrixStats::rowVars(sce_means),
                              row.names(sce_means))
    #### Convert to deciles ####
    deciles <- calc_quantiles(v = gene_variance,
                              n_quantiles = n_quantiles,
                              report_filters = verbose)
    #### Remove genes below the min_variance_decile ####
    gene_variance <- gene_variance[deciles>=min_variance_decile]
    # DelayedMatrixStats::rowQuantiles(gene_variance)
    # hist(log(gene_variance), breaks=50)
    printer(paste(nrow(sce)-length(gene_variance),"/",nrow(sce),
                  "genes dropped @ DGE variance_decile ≥",min_variance_decile), v=verbose)
    sce <- sce[names(gene_variance),]
    return(sce)
}



DelayedArray_grouped_stats <- function(sce,
                                       grouping_var,
                                       return_sce=T,
                                       stat="mean",
                                       rownames_var="Gene"){
    cdata <- SummarizedExperiment::colData(sce)
    if(!grouping_var %in% colnames(cdata)) stop(grouping_var," is not a column in colData(sce).")
    core_allocation <- assign_cores(worker_cores=.90, verbose=F)

    lv2_groups <- unique(cdata[[grouping_var]])
    sce <- check_sce_colnames(sce, colnames_var = grouping_var)
    sce <- check_sce_rownames(sce, rownames_var = rownames_var)

    exp_agg <- parallel::mclapply(lv2_groups, function(x){
        message_parallel("Calculating ",stat," : ",x)
        if(stat=="mean"){
            agg <- DelayedMatrixStats::rowMeans2(SummarizedExperiment::assay(sce),
                                                 cols = x, na.rm = T)
        }
        if(stat=="var"){
            agg <- DelayedMatrixStats::rowVars(SummarizedExperiment::assay(sce),
                                               cols = x, na.rm = T)
        }
        return(agg)
        # IMPORTANT!: don't parallelize here bc it causes conflicts with DelayedArray functions (when they use multiple chunks)
    }, mc.cores = 1 #core_allocation$worker_cores
    ) %>%
        `names<-`(lv2_groups) %>%
        data.table::as.data.table() %>%
        `rownames<-`(row.names(SummarizedExperiment::assay(sce))) %>%
        DelayedArray::DelayedArray()
    if(return_sce) {
        sce_agg <- SingleCellExperiment::SingleCellExperiment(
            assays      = list(raw = DelayedArray::DelayedArray(as(exp_agg,"sparseMatrix")) ),
            # colData     = S4Vectors::DataFrame(object$col.attrs[[1]]),
            rowData     =  SummarizedExperiment::rowData(sce)
        )
        # Slightly modify colData
        cdata <- SummarizedExperiment::colData(sce_agg)
        cdata[[grouping_var]] <- row.names(cdata)
        SummarizedExperiment::colData(sce_agg) <- cdata
        return(sce_agg)
    } else{ return(exp_agg) }
}


DelayedArray_normalize <- function(sce,
                                   log_norm=T,
                                   min_max=T,
                                   plot_hists=F){
    core_allocation <- assign_cores(worker_cores = .90)
    mat <- SummarizedExperiment::assay(sce)
    if(log_norm){
        mat_log <- log1p(mat)
        mat <- mat_log
    }
    if(min_max){
        col_max <-DelayedArray::colMaxs(mat, na.rm = T)
        col_min <-  DelayedArray::colMins(mat, na.rm = T)
        mat_normed <- DelayedArray::t( (DelayedArray::t(mat)- col_min) / (col_max - col_min) )
        mat <- mat_normed
    }
    if(plot_hists){
        hist(DelayedArray::colMeans(mat, na.rm = T)) %>% print()
        hist(DelayedArray::colMeans(mat_log, na.rm = T)) %>% print()
        hist(DelayedArray::colMeans(mat_normed, na.rm = T)) %>% print()
    }
    return(mat)
}



run_MAST <- function(sce){
    sca_delay <- MAST::FromMatrix(
        exprsArray = list(SummarizedExperiment::assay(sce)),
        cData = SummarizedExperiment::colData(sce),
        fData = SummarizedExperiment::rowData(sce),
        check_sanity = F
    )
    options(mc.cores=core_allocation$worker_cores)
    res =  MAST::zlm( ~ Class,
                      sca = sca_delay,
                      parallel = T)
    MAST::waldTest(res, MAST::CoefficientHypothesis('Stim.ConditionUnstim'))
}


run_limma <- function(sce,
                      level2annot,
                      verbose=T,
                      ...){
    messager("DGE:: Limma...",v=verbose)
    exp <- SummarizedExperiment::assay(sce)
    cdata <- SummarizedExperiment::colData(sce)
    ## Prepare groupings
    if(length(level2annot)==1){
        level2_options <- cdata[[level2annot]]
    } else { level2_options <- level2annot }
    level2_options <- as.factor(as.character(level2_options))

    mod_matrix  = model.matrix(~level2_options)
    fit = limma::lmFit(exp, mod_matrix, ...)
    eb = limma::eBayes(fit)
    return(eb)
}



run_DESeq2 <- function(sce,
                       level2annot,
                       verbose=T,
                       ...){
    messager("DGE:: DESeq2...",v=verbose)
    if(!"DESeq2" %in% row.names(installed.packages())){
        stop("Please install DESeq2 first: BiocManager::install('DESeq2')")
    }
    core_allocation <- assign_cores(worker_cores=.90)
    exp <- SummarizedExperiment::assay(sce)
    cdata <- SummarizedExperiment::colData(sce)
    # NOTE:: When you're running DESeq2 on sparse SCE data,
    ## there are two ways to avoid issues when DESeq() tries to log your data.
    ## 1) add 1 to you expression matrix (much faster).
    ## 2) set sfType = "iterate" to enable iterative size factor estimation (veerrrry slow).
    dds <- DESeq2::DESeqDataSetFromMatrix(exp+1,
                                          colData = cdata,
                                          design = formula(paste("~",level2annot)))
    dds <- DESeq2::DESeq(dds,
                         # Best for scRNAseq data.
                         test="LRT",
                         reduced = ~1,
                         # DESeq2 v1.31.10 (not yet released on BioC)
                         # now has glmGamPoi directly integrated directly!
                         ## https://github.com/mikelove/DESeq2/issues/29
                         ## default="parametric"
                         # fitType="glmGamPoi",
                         # sfType = "iterate",
                         parallel = T,
                         ...)
    dds_res <- DESeq2::results(dds)
    return(dds_res)
}



#### OG drop.uninformative.genes ####
# drop.uninformative.genes <- function(exp,level2annot){
#     if(class(exp[1,1])=="character"){
#         exp = as.matrix(exp)
#         storage.mode(exp) <- "numeric"
#     }
#     level2annot = as.character(level2annot)
#     summed = apply(exp,1,sum)
#     exp = exp[summed!=0,]
#     mod  = model.matrix(~level2annot)
#     fit = lmFit(exp,mod)
#     eb = eBayes(fit)
#     pF = p.adjust(eb$F.p.value,method="BH")
#     exp = exp[pF<0.00001,]
#     return(exp)
# }



