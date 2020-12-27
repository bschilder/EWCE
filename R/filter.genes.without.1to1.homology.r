#' filter.genes.without.1to1.homolog
#'
#' \code{filter.genes.without.1to1.homolog} Takes the filenames of celltype_data files, loads them, and drops any genes
#' which don't have a 1:1 homolog based on biomart. The new files are saved to disc, appending '_1to1only' to the file
#' tag of the previous file.
#'
#' @param filenames Array of filenames for sct_data saved as .Rda files
#' @return Array of filenames included the ones with only 1:1 homologs
#' @examples
#' # Load the single cell data
#' data(cortex_mrna)
#' expData = cortex_mrna$exp
#' expData = expData[1:500,] # Use only a subset to keep the example quick
#' l1=cortex_mrna$annot$level1class
#' l2=cortex_mrna$annot$level2class
#' annotLevels = list(level1class=l1,level2class=l2)
#' fNames_ALLCELLS = generate.celltype.data(exp=expData,annotLevels=annotLevels,groupName="allKImouse")
#' fNames_ALLCELLS = filter.genes.without.1to1.homolog(fNames_ALLCELLS)
#' @export
#' @import utils
#' @import stringr
filter.genes.without.1to1.homolog <- function(filenames){
    newFilenames = filenames
    #data("mouse_to_human_homologs")
    mouse_to_human_homologs <- EWCE::mouse_to_human_homologs
    orthologs = mouse_to_human_homologs
    mgi_1to1   = orthologs$MGI.symbol
    hgnc_1to1   = orthologs$HGNC.symbol
    for(ff in filenames){
        ctd <- loadRData(ff)
        sct_genes = rownames(ctd[[1]]$mean_exp)
        # If it's a m
        if(sum(sct_genes %in% orthologs$MGI.symbol) > sum(sct_genes %in% orthologs$HGNC.symbol)){
            symbol_1to1_in_sct = mgi_1to1[mgi_1to1 %in% sct_genes]
        }else{ symbol_1to1_in_sct = hgnc_1to1[hgnc_1to1 %in% sct_genes] }
        ctd[[1]]$mean_exp = ctd[[1]]$mean_exp[symbol_1to1_in_sct,]
        ctd[[1]]$specificity = ctd[[1]]$specificity[symbol_1to1_in_sct,]
        #ff2 = gsub("_level","_1to1only_level",ff)
        ff2 = gsub("\\.rda","_1to1only\\.rda",ff)
        save(ctd,file=ff2)
        newFilenames = c(newFilenames,ff2)
    }
    return(newFilenames)
}






#' filter_nonorthologues
#'
#' \code{filter_nonorthologues} Takes the filenames of CellTypeData files, loads them, and drops any genes
#' which don't have a 1:1 orthologues based on \pkg{homologene}.
#' The new files are saved to disk, appending '_orthologues' to the file name.
#'
#' @param filenames List of file names for sct_data saved as \emph{.rda} files.
#' @param suffix Suffix to add to the file name (right before \emph{.rda}).
#' @inheritParams convert_orthologues
#' @return List of the filtered CellTypeData file names.
#' @examples
#' # Load the single cell data
#' data(cortex_mrna)
#' expData = cortex_mrna$exp
#' expData = expData[1:500,] # Use only a subset to keep the example quick
#' l1=cortex_mrna$annot$level1class
#' l2=cortex_mrna$annot$level2class
#' annotLevels = list(level1class=l1,level2class=l2)
#' fNames_ALLCELLS = generate.celltype.data(exp=expData,annotLevels=annotLevels,groupName="allKImouse")
#' fNames_ALLCELLS_orths = filter_nonorthologues(fNames_ALLCELLS)
#' @export
#' @import utils
#' @import stringr
filter_nonorthologues <- function(filenames,
                                  input_species="mouse",
                                  annot_levels=NULL,
                                  suffix="_orthologues",
                                  verbose=T){
    newFilenames <- filenames
    for(ff in filenames){
        ctd <- loadRData(ff)
        annot_levels <- if(is.null(annot_levels)) 1:length(ctd)
        for(lvl in annot_levels){
            messager("+ Processing level",lvl,"...", v=verbose)
            orths <- convert_orthologues(gene_df=ctd[[lvl]]$mean_exp,
                                         gene_col="rownames",
                                         input_species=input_species,
                                         drop_nonhuman_genes=T,
                                         one_to_one_only=T,
                                         genes_as_rownames=T,
                                         verbose=verbose)
            ctd[[lvl]]$mean_exp = ctd[[lvl]]$mean_exp[orths$Gene_orig,]
            if("specificity" %in% names(ctd[[lvl]])){
                ctd[[lvl]]$specificity = ctd[[lvl]]$specificity[orths$Gene_orig,]
            }
            if("specificity_quantiles" %in% names(ctd[[lvl]]) ){
                ctd[[lvl]]$specificity_quantiles = ctd[[lvl]]$specificity_quantiles[orths$Gene_orig,]
            }
        }
        ff2 = gsub("\\.rda",paste0(suffix,"\\.rda"),ff)
        save(ctd, file=ff2)
        newFilenames = c(newFilenames,ff2)
    }
    return(newFilenames)
}











#' Convert genes from one species to another
#'
#' Currently supports: mouse ==> human, fly ==> human.
#' @param gene_df Data.frame containing the genes symbols of the input species.
#' @param  gene_col Character string indicating either the input species' gene symbols are in a
#' column with the input species gene symbols (\code{gene_col="<gene_col_name>"})
#' or are the rownames (\code{gene_col="rownames"}).
#' @param input_species Name of the input species ("mouse" or "fly").
#' @param drop_nonhuman_genes Drop genes that don't have an orthologue in humans.
#' @param one_to_one_only Drop genes that don't have a 1:1 mapping between input species and human
#' (i.e. drop genes with multiple human orthologues).
#' @param genes_as_rownames Whether to return the data.frame with the row names set to the human orthologues.
#' @param verbose Show progress messages.
#' @return Data.frame with the gene symbols of the input species ("Gene_orig")
#' and converted human orthologues ("Gene").
#' @import homologene
convert_orthologues <- function(gene_df,
                                gene_col="Gene",
                                input_species="mouse",
                                drop_nonhuman_genes=F,
                                one_to_one_only=T,
                                genes_as_rownames=F,
                                verbose=T){
    gene_df <- data.frame(gene_df)
    messager("Converting genes: ",input_species," ===> human", v=verbose)

    if(gene_col %in% c("rownames","row.names")){
        printer("+ Converting rownames to Gene col...",v=verbose)
        gene_df$Gene <- row.names(gene_df)
        gene_col <- "Gene"
    }
    if("Gene_orig" %in% colnames(gene_df)){
        printer("+ Removing Gene_orig col...",v=verbose)
        gene_df <- dplyr::select(gene_df, -Gene_orig)
    }

    printer("+ Searching for orthologues...",v=verbose)
    taxaID_dict <- function(species=c("human","mouse","fly")){
        dict <- c("human"=9606,
                  "mouse"=10090,
                  "fly"=7227)
        return(dict[species])
    }
    taxaID <- taxaID_dict(species=input_species)
    input_genes <- gene_df[[gene_col]]
    orths <- homologene::homologene(genes = input_genes,
                                    inTax = unname(taxaID),
                                    outTax = 9606)
    orths <- orths[orths[,1] %in% input_genes,]
    orths_key <- setNames(orths[,2], orths[,1])
    gene_df <- dplyr::rename(gene_df, Gene_orig=all_of(gene_col))
    gene_df$Gene <- orths_key[input_genes]

    if(drop_nonhuman_genes){
        printer("+ Dropping genes with no human orthologues...",v=verbose)
        gene_df <- gene_df[!is.na(gene_df$Gene),]
    }

    if(one_to_one_only){
        printer("+ Dropping genes that don't have 1-to-1 gene mappings...",v=verbose)
        gene_df <- gene_df[!duplicated(gene_df$Gene_orig),]
        gene_df <- gene_df[!duplicated(gene_df$Gene),]
    }
    if(genes_as_rownames){
        if(drop_nonhuman_genes==F){
            warning("+ drop_nonhuman_genes must =T in order to set genes_as_rownames=T.")
        }else {
            printer("+ Setting Gene col as rownames...",v=verbose)
            if(one_to_one_only==F) stop("Genes must be unique to be row.names. Try again with`one_to_one_only=T`")
            row.names(gene_df) <- gene_df$Gene
        }
    }
    # Report
    n_dropped <- length(unique(input_genes))-length(unique(gene_df$Gene))
    messager("Genes dropped during inter-species conversion: ",
            format(n_dropped,big.mark=","),
            " / ",
            format(length(unique(input_genes)),big.mark=","),
            " (",format(n_dropped/length(unique(input_genes))*100, digits = 4),"%)",
            v=verbose)
    return(gene_df)
}

