#' generate.celltype.data
#'
#' \code{generate.celltype.data} Takes expression & cell type annotations and creates celltype_data files which contain the mean and specificity matrices
#'
#' @param exp Numerical matrix with row for each gene and column for each cell. Row names are MGI/HGNC gene symbols.
#' Column names are cell IDs which can be cross referenced against the annot data frame.
#' @param annotLevels List with arrays of strings containing the cell type names associated with each column in exp.
#' Alternatively, if \code{exp} is of class \code{SingleCellExperiment}, \code{annotLevels}.
#' can simply be a list containing the names of the columns you wish to use as levels (so long as they are in in \code{colnames(exp@@colData))}).
#' @param groupName A human readable name for referring to the dataset being loaded.
#' @param no_cores Number of cores that should be used to speedup the computation. Use no_cores = 1 when using this package in windows system.
#' @param savePath Directory where the CTD file should be saved.
#' @param force_new_file If a file of the same name as the one being created already exists, overwrite it.
#' @param return_ctd Return the CTD object in a list along with the file name, instead of just the file name.
#' @return File names for the saved CellTypeData (CTD) files.
#' @examples
#' # Load the single cell data
#' data("cortex_mrna")
#'
#' #### Run as Sparse Matrix (original method) ####
#' expData = cortex_mrna$exp
#' expData = expData[1:500,] # Use only a subset to keep the example quick
#' l1=cortex_mrna$annot$level1class
#' l2=cortex_mrna$annot$level2class
#' annotLevels = list(l1=l1,l2=l2)
#' fNames_ALLCELLS = generate.celltype.data(exp=expData,annotLevels,"allKImouse")
#'
#'
#' #### Run as DelayedArray (for large datasets) ####
#' exp <- cortex_mrna$exp
#' colData <- cortex_mrna$annot
#' save_dir <- "./cortex_mrna_HDF5Array"
#' sce <- convert_to_SCE(exp=exp, colData=colData, save_dir=save_dir)
#' print(sce)
#'
#' annotLevels <- c("level1class","level2class")
#' fNames_SCE <-generate.celltype.data(exp=sce, annotLevels=annotLevels, groupName="mouse_cortex_SCE", savePath="./")
#' @export
#' @import parallel
#' @import future
#' @import ggdendro
#' @import gridExtra
#' @importFrom Matrix Matrix
#' @import RNOmni
#' @import ggdendro
#' @source
#' #' \href{https://petehaitch.github.io/BioC2020_DelayedArray_workshop/articles/Effectively_using_the_DelayedArray_framework_for_users.html}{DelayedArray workshop}
generate.celltype.data <- function(exp,
                                   annotLevels,
                                   groupName,
                                   numberOfBins=40,
                                   no_cores=NULL,
                                   add_names=F,
                                   savePath="~/",
                                   file_prefix="CellTypeData",
                                   force_new_file=T,
                                   return_ctd=F,
                                   drop_nonhuman_genes=F,
                                   input_species="mouse",
                                   verbose=T){
    #### Check group name ####
    if(is.null(groupName)){stop("ERROR: groupName must be set. groupName is used to label the files created by this function.")}
    if(groupName==""){stop("ERROR: groupName must be set. groupName is used to label the files created by this function.")}

    #### Check if file already exists ####
    fNames <- sprintf("%s/%s_%s.rda",savePath,file_prefix,groupName)
    if(file.exists(fNames) & force_new_file==F){
        messager("+ Pre-existing file of the same name detected. Use `force_new_file=T` to overwrite this file.",v=verbose)
        messager("+ Returning pre-existing file path.",v=verbose)
        return(fNames)
    }

    #### Prepare cores ####
    worker_cores <- if(is.null(no_cores)) .90 else no_cores
    core_allocation <- assign_cores(worker_cores = worker_cores)
    no_cores <- core_allocation$worker_cores

    # Calculate summary stats from matrix

    #### DelayedArray method ####
    if(class(exp)[1] %in% c("SingleCellExperiment","SingleCellLoomExperiment")){
        messager("Processing as DelayedArray...",v=verbose)
        sce <- exp

        if(drop_nonhuman_genes){
            sce <- check_sce_rownames(sce, rownames_var = "Gene")
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
        #### Prepare annotation levels ####
        if(length(annotLevels[[1]])>1){
            level_names <- names(annotLevels)
        }else {
            level_names <- annotLevels
        }
        if(any(!level_names %in% colnames(SummarizedExperiment::colData(sce))) ){
            stop("+ When exp is of class 'SingleCellExperiment'",
                 " all names(annotLevels) must be in the column names of exp@colData.")
        }

        # IMPORTANT!: Do NOT parallelize this step. It will interfere with DelayedArray functions.
        ctd <- lapply(level_names,function(lvl, .sce=sce){
            messager("+ Processsing level = ",lvl,v=verbose)
            grouped_means <- DelayedArray_grouped_stats(sce = .sce,
                                                        grouping_var = lvl,
                                                        stat = "mean",
                                                        return_sce = F)
            grouped_specificity <- grouped_means / DelayedArray::rowSums(grouped_means)
            # Convert to sparse DelayedMatrices
            grouped_means <- DelayedArray::DelayedArray(as(grouped_means, "sparseMatrix"))
            grouped_specificity <- DelayedArray::DelayedArray(as(grouped_specificity, "sparseMatrix"))
            return(list(mean_exp=grouped_means,
                        specificity=grouped_specificity))
        })
        if(add_names) names(ctd) <- level_names
        # Convert the DelayedMatrices as sparseMatrices
        # ctd_sparse <- ctd
        # for(lvl in names(ctd_sparse)){
        #     printer(lvl)
        #     ctd_lvl <- ctd_sparse[[lvl]]
        #     for(metric in names(ctd_lvl)){
        #         printer("+",metric)
        #         ctd_lvl[[metric]] <- DelayedArray::DelayedArray(as(ctd_lvl[[metric]], "sparseMatrix"))
        #     }
        #     ctd_sparse[[lvl]] <- ctd_lvl
        # }
        # saveRDS(ctd_sparse,"../model_celltype_conservation/processed_data/EWCE/LaManno2020/ctd_LaManno2020_4levels_sparseDA.RDS")


        # ***** DONE ***** #
    } else {
        #### Original method ####
        messager("Processing in-memory...",v=verbose)

        if(sum(is.na(exp))>0){stop("NA values detected in expresson matrix. All NA values should be removed before calling EWCE.")}
        # First, check the number of annotations equals the number of columns in the expression data
        out <- lapply(annotLevels,test <- function(x,exp){if(length(x)!=dim(exp)[2]){stop("Error: length of all annotation levels must equal the number of columns in exp matrix")}},exp)

        # Initialize CTD list
        ctd = list()
        for(i in 1:length(annotLevels)){ctd[[length(ctd)+1]] = list(annot=annotLevels[[i]])}

        # Convert characters to numbers
        matrix_classes <- c("data.table","data.frame","tbl_df","tbl","matrix","Matrix","array","DelayedArray","DelayedMatrix",names(getClass("Matrix")@subclasses)) # more than 40 ..
        if(!class(exp)[1] %in% matrix_classes){
            printer("Converting exp to matrix (double)",v=verbose)
            exp <- suppressWarnings(apply(exp,2,function(x) {storage.mode(x) <- 'double'; x}))
        }

        # Make exp into a sparse matrix
        ## Matrix() will make data sparse automatically if >50% of the data is 0s.
        ## But we set sparse=T here just to be explicit.
        printer("Converting exp to sparse DelayedMatrix",v=verbose)
        exp <- DelayedArray::DelayedArray(as(exp, "sparseMatrix"))

        calculate.meanexp.for.level <- function(ctd_oneLevel,
                                                expMatrix,
                                                verbose=T){
            if( dim(expMatrix)[2]==length(unique(ctd_oneLevel$annot)) ){
                print(dim(expMatrix)[2])
                print(length(ctd_oneLevel$annot))
                if(sum(!colnames(expMatrix)==ctd_oneLevel$annot)!=0){
                    stop("There are an equal number of celltypes in expMatrix and ctd_oneLevel but the names do not match")
                }
                messager("+ Assuming supplied matrix is already a mean expression matrix. Skipping calculation.",v=verbose)
                ctd_oneLevel$mean_exp = as.matrix(expMatrix)
            }else{
                # Sum reads in each cell type
                #mean_exp = apply(expMatrix,1,aggregate.over.celltypes,ctd_oneLevel$annot)
                mm <- model.matrix(~ 0 + ctd_oneLevel$annot)
                colnames(mm) <- names(table(ctd_oneLevel$annot))
                mat.summary.mm1 <- expMatrix %*% mm

                # Divide by the number of cells to get the mean
                cellCounts = table(ctd_oneLevel$annot)
                for(i in 1:dim(mat.summary.mm1)[2]){mat.summary.mm1[,i] = mat.summary.mm1[,i]/cellCounts[i]}

                ctd_oneLevel$mean_exp = as.matrix(mat.summary.mm1)
            }
            return(ctd_oneLevel)
        }

        calculate.specificity.for.level <- function(ctd_oneLevel){
            normalised_meanExp = t(t(ctd_oneLevel$mean_exp)*(1/colSums(ctd_oneLevel$mean_exp)))
            ctd_oneLevel$specificity = normalised_meanExp/(apply(normalised_meanExp,1,sum)+0.000000000001)
            return(ctd_oneLevel)
        }
        printer("+ Calculating normalized mean expression.",v=verbose)
        ctd2 = parallel::mclapply(ctd,calculate.meanexp.for.level,exp,verbose,mc.cores=no_cores)
        printer("+ Calculating normalized specificity.",v=verbose)
        ctd3 = parallel::mclapply(ctd2,calculate.specificity.for.level,mc.cores=no_cores)
        ctd=ctd3
        #stopCluster(cl)
        # ***** DONE ***** #
    }

    #### Quantile normalization ####
    # Use the rank norm transformation on specificity
    # rNorm <- function(ctdIN){   bbb = t(apply(ctdIN$specificity,1,RNOmni::rankNorm));  return(bbb)    }

    # ADD DENDROGRAM DATA TO CTD
    ctd = lapply(ctd, bin.specificity.into.quantiles, numberOfBins=numberOfBins)
    ctd = lapply(ctd, prep.dendro)

    #### Save results ####
    messager("+ Saving results ==> ",fNames, v=verbose)
    dir.create(dirname(fNames), showWarnings = F, recursive = T)
    save(ctd, file=fNames)

    if(return_ctd){
        messager("+ Returning list of CTD file name, and the CTD itself.",v=verbose)
        return(list(fNames=fNames,
                    ctd=ctd))
    } else { return(fNames)}
}







grouped_col_means <- function(sce,
                              group_var="Class",
                              verbose=T){
    # library(DelayedArray); library(dplyr);
    DelayedArray:::set_verbose_block_processing(verbose)
    grouped_sums <- DelayedArray::colsum(sce@assays@data$raw,
                                         group = sce[[group_var]],
                                         reorder = F) %>%
        `rownames<-`(row.names(sce))
    grouped_colcounts <- as.vector(table(sce[[group_var]]))
    grouped_means <- grouped_sums / grouped_colcounts
    return(grouped_means)
}

grouped_col_specificity <- function(sce,
                                    group_var="Class",
                                    verbose=T){
    # library(DelayedArray)
    grouped_means <- grouped_col_means(sce,
                                       group_var=group_var,
                                       verbose=T)
    grouped_specificity <- grouped_means / rowSums(grouped_means)
    return(grouped_specificity)
}
