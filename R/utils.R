
#' Make pretty tables
#'
#' Take any data.frame and turn it into an interactive, scrollable table
#' with buttons to download the table in multiple formats.
#'
#' @export
create_DT <- function (DF, caption = "", scrollY = 400)
{
    data <- DT::datatable(DF, caption = caption, extensions = "Buttons",
                          options = list(dom = "Bfrtip", buttons = c("copy", "csv",
                                                                     "excel", "pdf", "print"), scrollY = scrollY, scrollX = T,
                                         scrollCollapse = T, paging = F, columnDefs = list(list(className = "dt-center",
                                                                                                targets = "_all"))))
    return(data)
}




message_parallel <- function(...){
    system(sprintf('echo "%s"', paste0(..., collapse="")))
}

messager <- function(..., v=T){if(v){message(paste(...))}}

printer <- function(..., v=T){if(v){print(paste(...))}}


#' loadRData
#'
#' Load processed data using a function that assigns it
#' to a specific variable (so you don't have to guess what the loaded variable name is).
#'
#' @export
loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
}




assign_cores <- function(worker_cores=.90,
                         verbose=T){
    # Enable parallelization of HDF5 functions
    ## Allocate ~10% of your available cores to non-parallelized processes
    total_cores <- parallel::detectCores()
    if(worker_cores<1){
        reserved_cores <-  ceiling(total_cores*(1-worker_cores))
        workers <- total_cores - reserved_cores
    } else {
        workers <- worker_cores
        reserved_cores <-  total_cores - workers
    }
    message("+ ",workers," core(s) assigned as workers (",reserved_cores," reserved).")
    DelayedArray::setAutoBPPARAM(BiocParallel::MulticoreParam(workers))
    DelayedArray:::set_verbose_block_processing(verbose)
    return(list(worker_cores=workers,
                reserved_cores=reserved_cores,
                total_cores=total_cores))
}



calc_quantiles <- function(v,
                           n_quantiles=10,
                           report_filters=T){
    calc_percentile <- stats::ecdf(1:n_quantiles)
    quantiles <- calc_percentile(v)
    # Report the number of items left after filtering at each quantile
    if(report_filters){
        print("N remaining at each quantile filter:")
        genes_left <- lapply(unique(quantiles), function(x){
            setNames(length(v[quantiles>=x]), paste0("≥",x))
        }) %>% unlist() %>% sort() %>% rev()
        print(genes_left)
    }
    return(quantiles)
}

sce_filepath <- function(sce,
                         sce_save_dir=NULL){
    file_info <- DelayedArray::seed(SummarizedExperiment::assay(sce))
    if("filepath" %in% slotNames(file_info)){
        return(file_info@filepath)
    } else { return(NULL) }
}
