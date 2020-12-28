#' drop.uninformative.genes
#'
#' \code{drop.uninformative.genes} drops genes from an SCT expression matrix if they do not significantly vary between any celltypes.
#' Makes this decision based on use of an ANOVA (implemented with Limma). If the F-statistic for variation amongst type2 annotations
#' is less than a strict p-threshold, then the gene is dropped.
#'
#' @param exp Expression matrix with gene names as rownames.
#' @param level2annot Array of cell types, with each sequentially corresponding a column in the expression matrix
#' @return exp Expression matrix with gene names as rownames.
#' @examples
#' data("cortex_mrna")
#' cortex_mrna$exp = cortex_mrna$exp[1:300,] # Use only a subset of genes to keep the example quick
#' exp2 = drop.uninformative.genes(exp=cortex_mrna$exp, level2annot=cortex_mrna$annot$level2class)
#' @export
#' @import limma
#' @import stats
drop.uninformative.genes <- function(exp,
                                     level2annot,
                                     as_DelayedArray=F,
                                     colData=NULL,
                                     rowData=NULL,
                                     sce_save_dir=NULL,
                                     pseudobulk_by=NULL,
                                     adj_pval_thresh=0.00001,
                                     verbose=T){
    # Allow user to supply either the vector
    ## or simply the name of the column in the colData
    if(as_DelayedArray){
        #### DelayedArray method ####
        # Best for large datasets that can't fit into memory.
        message("Processing as DelayedArray...")

        if(class(exp)[1]!="SingleCellExperiment"){
            sce <- construct_SCE(exp=exp,
                                  colData=colData,
                                  rowData=rowData,
                                  save_dir=sce_save_dir,
                                  verbose=verbose)
        }else {message("+ `exp` is of class 'SingleCellExperiment'"); sce <- exp}

        if((length(level2annot)>1) | (!level2annot %in% colnames(sce@colData)) ){
            # Required because `construct_SCE()` does some filtering of cells
            ## (thus the original user-supplied vector (as level2annot) might not line up with the exp data any more).
            stop("+ `level2annot` must be the name of a column in colData (e.g. level2annot='celltype').")
        }
        message("+ Running Gamma-Poisson Generalized Linear regression (glmGamPoi)...")
        sce_de <- run_glmGamPoi_DE(sce,
                                   level2annot,
                                   pseudobulk_by=pseudobulk_by,
                                   pval_adjust_method="BH",
                                   adj_pval_thresh=adj_pval_thresh,
                                   on_disk=T,
                                   return_as_SCE=T,
                                   verbose=verbose)
        return(sce_de)
    } else {
        #### Original method ####
        # Best for smaller datasets that can fit into memory)
        message("Processing in-memory...")
        if(class(exp[1,1])=="character"){
            exp = as.matrix(exp)
            storage.mode(exp) <- "numeric"
        }
        if(length(level2annot)==1 & (!is.null(colData)) ){
            level2annot <- colData[[level2annot]]
        }
        level2annot = as.factor(as.character(level2annot))
        summed = apply(exp,1,sum)
        exp = exp[summed!=0,]
        mod_matrix  = model.matrix(~level2annot)
        fit = limma::lmFit(exp, mod_matrix)
        eb = limma::eBayes(fit)
        pF = stats::p.adjust(eb$F.p.value,method="BH")
        exp = exp[pF<0.00001,]
        return(exp)
    }
}




#' Construct SingleCellExperiment
#'
#' @import DelayedArray
#' @import BiocParallel
#' @import HDF5Array
#' @examples
#' data("cortex_mrna")
#' exp <- cortex_mrna$exp
#' colData <- cortex_mrna$annot
#' save_dir <- "./cortex_mrna_HDF5Array"
#' sce <- construct_SCE(exp=exp, colData=colData, save_dir=save_dir)
#' print(sce)
#'  @export
#'  @import DelayedArray
#'  @import BiocParallel
#'  @import parallel
#'  @import SingleCellExperiment
construct_SCE <- function(exp,
                           colData=NULL,
                           rowData=NULL,
                           save_dir=NULL,
                           drop_empty=T,
                           replace_HDF5=F,
                           quicksave_HDF5=F,
                           verbose=T){
    # library(DelayedArray)
    # library(BiocParallel)
    # library(SingleCellExperiment)
    core_allocation <- assign_cores(worker_cores = .90)

    message("+ Constructing SingleCellExperiment...")
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays      = list(raw = DelayedArray::DelayedArray(exp)),
        colData     = colData,
        rowData     = rowData
    )
    if(drop_empty){
        message("+ Dropping rows and/or cols with sums==0...")
        non_empty_rows <- which(DelayedMatrixStats::rowSums2(SummarizedExperiment::assay(sce)) > 0)
        non_empty_cols <- which(DelayedMatrixStats::colSums2(SummarizedExperiment::assay(sce)) > 0)
        # sce_sub <- sce[sample(non_empty_rows,300), non_empty_cols]
        sce <- sce[non_empty_rows, non_empty_cols]
    }
    sce <- save_SCE(sce=sce,
                    save_dir=save_dir,
                    quicksave_HDF5=quicksave_HDF5,
                    replace_HDF5=replace_HDF5,
                    verbose=verbose)
    return(sce)
}





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
                             verbose=T){
    level2annot <- as.factor(sce[[level2annot]])
    mod_matrix  <- model.matrix(~level2annot)
    fit <- glmGamPoi::glm_gp(sce,
                             design = mod_matrix,
                             # on_disk=F when testing on subsets of SCE
                             on_disk = on_disk,
                             verbose = verbose)
    # Save fitted model as intermediate
    sce_dir <- dirname(DelayedArray::seed(SummarizedExperiment::assay(sce))@filepath)
    fit_path <- file.path(sce_dir,paste(basename(sce_dir),"glm_gp.RDS",sep="."))
    messager("+ Saving intermediate file ==>",fit_path)
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
    message("+ ",genes_dropped," / ",nrow(sce),
            " (",round(genes_dropped/nrow(sce)*100, 1),"%)",
            " genes dropped due to lack of differential expression ",
            "between level2 annotations.")
    if(return_as_SCE){
        return(sce_de)
    }else {
        return(as(SummarizedExperiment::assay(sce_de), "matrix"))
    }
}



