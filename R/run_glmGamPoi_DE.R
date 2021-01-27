
#' Run Gamma-Poisson Generalized Linear regression (glmGamPoi)
#'
#' Installation troubleshooting:
#' 1. Try to install using  Bioconductor: \code{BiocManager::install("glmGamPoi")}
#' 2. If that fails, try installing directly from GitHub: \code{devtools::install_github("const-ae/glmGamPoi")}
#' 3. If that fails due to an error like \code{ld: warning: directory not found for option '-L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin15/6.1.0'}
#' follow \href{https://stackoverflow.com/questions/35999874/mac-os-x-r-error-ld-warning-directory-not-found-for-option}{these instructions} to fix this.
#' @examples
#' data("cortex_mrna")
#' exp <- cortex_mrna$exp[1:300,]
#' colData <- cortex_mrna$annot
#' save_dir <- "~/Desktop/cortex_mrna_HDF5Array"
#' sce <- construct_SCE(exp=exp, colData=colData, save_dir=save_dir, replace_HDF5 = T)
#' exp <- run_glmGamPoi_DE(sce, level2annot="level2class")
#' @import glmGamPoi
run_glmGamPoi_DE <- function(sce,
                             level2annot,
                             pseudobulk_by=NULL,
                             pval_adjust_method="BH",
                             adj_pval_thresh=0.00001,
                             on_disk=T,
                             update_on_disk=T,
                             return_as_SCE=T,
                             sce_save_dir=NULL,
                             verbose=T,
                             ...){
    core_allocation <- assign_cores(worker_cores=.90, verbose=F)
    # You must make sure the sce is
    # saved to disk in order to get its filepath.
    h5_path <- sce_filepath(sce)
    if(is.null(h5_path)){
        if(is.null(sce_save_dir)){
            stop("`sce_save_dir` required to save SCE object to disk.")
        }else {
            sce <- HDF5Array::saveHDF5SummarizedExperiment(sce,
                                                           dir=sce_save_dir,
                                                           replace = T)
            h5_path <- sce_filepath(sce)
        }
    }


    try({
        # Make sure on_disk version matches the sce pointer
        sce <- HDF5Array::saveHDF5SummarizedExperiment(sce,
                                                       dir=sce_save_dir,
                                                       replace = T)
    })
    try({
        sce <- HDF5Array::quickResaveHDF5SummarizedExperiment(sce, verbose = verbose)
    })

    level2_options <- as.factor(sce[[level2annot]])
    mod_matrix  <- model.matrix(~level2_options)
    fit <- glmGamPoi::glm_gp(sce,
                             design = mod_matrix,
                             # on_disk=F when testing on subsets of SCE
                             on_disk = on_disk,
                             # Offset necessary for sparse scRNAseq data?
                             # offset = 1,
                             verbose = verbose)
    # Save fitted model as intermediate
    sce_dir <- dirname(h5_path)
    fit_path <- file.path(sce_dir,paste(basename(sce_dir),"glm_gp.RDS",sep="."))
    messager("+ Saving intermediate file ==>",fit_path, v=verbose)
    saveRDS(fit, fit_path)

    # Run DGE
    # intercept <- colnames(fit$Beta)[1]
    # normed <- (fit$Beta[,1] / sum(fit$Beta[,1]))
    # normed_scaled <- scales::rescale(normed, to = c(0,1))
    # fit$Beta <- cbind(fit$Beta, Intercept_normed=normed_scaled)
    de_res <- glmGamPoi::test_de(fit,
                                 # `pseudobulk_by` reduces false positives drastically
                                 pseudobulk_by = pseudobulk_by,
                                 contrast =`(Intercept)`,
                                 pval_adjust_method = pval_adjust_method,
                                 verbose = verbose)
    # Add DGE results back into SCE object
    sce_de <- SingleCellExperiment::SingleCellExperiment(
        assays      = list(raw = assay(sce)),
        colData     = sce@colData,
        rowData     = de_res
    )

    if(update_on_disk){
        messager("+ Updating the SCE object with the DGE results dataframe added to it...")
        sce_de <- HDF5Array::quickResaveHDF5SummarizedExperiment(sce_de, verbose=verbose)
    }

    # Only return deferentially expressed genes
    sce_de <- subset(sce_de, adj_pval<adj_pval_thresh)
    genes_dropped <- nrow(sce)-nrow(sce_de)
    messager(genes_dropped,"/",nrow(sce),
            "(",round(genes_dropped/nrow(sce)*100, 1),"%)",
            "genes dropped @ DGE adj_pval_thresh <",adj_pval_thresh, v=verbose)
    if(return_as_SCE){
        return(sce_de)
    }else {
        return(as(SummarizedExperiment::assay(sce_de), "matrix"))
    }
}
