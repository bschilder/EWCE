Expression Weighted Celltype Enrichment with *EWCE*
================
Nathan Skene & Brian Schilder
2021-02-11

<!-- Readme.md is generated from Readme.Rmd. Please edit that file -->

<!-- badges: start -->

[![Build
Status](https://travis-ci.org/NathanSkene/EWCE.svg?branch=master)](https://travis-ci.org/NathanSkene/EWCE)
<!-- badges: end -->

# Introduction

The *EWCE* package is designed to facilitate expression weighted
celltype enrichment analysis as described in our Frontiers in
Neuroscience paper<sup>1</sup>.

The package was originally designed to work with the single cell
cortical transcriptome data from the Linnarsson lab<sup>2</sup> which is
available at <http://linnarssonlab.org/cortex/>. Using this package it
is possible to read in any single cell transcriptome data, provided that
you have a cell by gene expression matrix (with each cell as a seperate
column) and a seperate annotation dataframe, with a row for each cell.

The *EWCE* process involves testing for whether the genes in a target
list have higher levels of expression in a given cell type than can
reasonably be expected by chance. The probability distribution for this
is estimated by randomly generating gene lists of equal length from a
set of background genes.

The *EWCE* method can be applied to any gene list. In the paper we
reported it’s application to genetic and transcriptomic datasets, and in
this vignette we detail how this can be done.

Note that throughout this vignette we use the terms ‘cell type’ and
‘sub-cell type’ to refer to two levels of annotation of what a cell
type is. This is described in further detail in our paper<sup>1</sup>,
but relates to the two levels of annotation provided in the Linnarsson
dataset<sup>2</sup>. In this dataset a cell is described as having a
cell type (i.e. ‘Interneuron’) and subcell type (i.e. ‘Int11’ a.k.a
Interneuron type 11).

# Overview

The process for using *EWCE* essentially involves three steps.

First, one needs to load the relevant single cell transcriptome dataset.
Single cell transcriptome data is read in from a text file using the
`read_celltype_data`.

The user then obtains a gene set and a suitable background gene set. As
the choice of gene sets is up to the user we do not provide functions
for doing this. Appropriate choice of background set is discussed in the
associated publication.

Bootstrapping is then performed using the `bootstrap.enrichment.test`
function.

# Installing EWCE

The *EWCE* package is available from github. The old version available
from the Bioconductor archives is depreciated and should not be used. To
be able to install the package one needs first to install R then run the
following lines of code:

    if (!require("devtools")) {
      install.packages("devtools")
    }
    devtools::install_github("nathanskene/ewce")

You can then load the package:

``` r
library(EWCE)
```

# Using with docker

Images with the latest version of EWCE are regularly pushed to
Dockerhub. If you already have docker installed you can load up a
working copy using the following commands. Note, that you will need to
replace the initial directory path with a location on your computer that
you wish to be able to access from within the docker image.

    docker pull nathanskene/ewce
    docker run --name=ewce -e PASSWORD=ewcedocker -p 8790:8790 -d -v /User/$USER:/var/ewce nathanskene/ewce:latest
    docker exec -ti ewce R

# Getting started

See the [vignette
website](https://nathanskene.github.io/EWCE/articles/EWCE.html) for
up-to-date instructions on usage.

If you have any problems please do file an issue here on github.

# Citation

If you use the EWCE package as well then please cite

[Skene, et al. Identification of Vulnerable Cell Types in Major Brain
Disorders Using Single Cell Transcriptomes and Expression Weighted Cell
Type Enrichment. Front. Neurosci,
2016.](https://www.frontiersin.org/articles/10.3389/fnins.2016.00016/full)

If you use the cortex/hippocampus single cell data associated with this
package then please cite the following papers:

[Zeisel, et al. Cell types in the mouse cortex and hippocampus revealed
by single-cell RNA-seq. Science,
2015.](http://www.sciencemag.org/content/early/2015/02/18/science.aaa1934.abstract)

# References

# Session Info

<details>

``` r
utils::sessionInfo()
```

    ## R version 4.0.2 (2020-06-22)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 20.10
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.10.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
    ##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
    ##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] EWCE_0.99.2
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] MatrixGenerics_1.2.1        Biobase_2.50.0             
    ##  [3] ggdendro_0.1.22             httr_1.4.2                 
    ##  [5] bit64_4.0.5                 assertthat_0.2.1           
    ##  [7] askpass_1.1                 stats4_4.0.2               
    ##  [9] BiocFileCache_1.14.0        blob_1.2.1                 
    ## [11] GenomeInfoDbData_1.2.4      yaml_2.2.1                 
    ## [13] progress_1.2.2              globals_0.14.0             
    ## [15] pillar_1.4.7                RSQLite_2.2.3              
    ## [17] lattice_0.20-41             limma_3.46.0               
    ## [19] glue_1.4.2                  digest_0.6.27              
    ## [21] XVector_0.30.0              GenomicRanges_1.42.0       
    ## [23] colorspace_2.0-0            htmltools_0.5.1.1          
    ## [25] Matrix_1.2-18               plyr_1.8.6                 
    ## [27] XML_3.99-0.5                pkgconfig_2.0.3            
    ## [29] biomaRt_2.46.3              listenv_0.8.0              
    ## [31] zlibbioc_1.36.0             purrr_0.3.4                
    ## [33] scales_1.1.1                tibble_3.0.6               
    ## [35] openssl_1.4.3               generics_0.1.0             
    ## [37] IRanges_2.24.1              ggplot2_3.3.3              
    ## [39] ellipsis_0.3.1              cachem_1.0.3               
    ## [41] SummarizedExperiment_1.20.0 BiocGenerics_0.36.0        
    ## [43] RNOmni_1.0.0                magrittr_2.0.1             
    ## [45] crayon_1.4.1                memoise_2.0.0              
    ## [47] evaluate_0.14               future_1.21.0              
    ## [49] parallelly_1.23.0           MASS_7.3-52                
    ## [51] homologene_1.4.68.19.3.27   xml2_1.3.2                 
    ## [53] tools_4.0.2                 prettyunits_1.1.1          
    ## [55] hms_1.0.0                   lifecycle_0.2.0            
    ## [57] matrixStats_0.58.0          stringr_1.4.0              
    ## [59] S4Vectors_0.28.1            munsell_0.5.0              
    ## [61] DelayedArray_0.16.1         AnnotationDbi_1.52.0       
    ## [63] compiler_4.0.2              GenomeInfoDb_1.26.2        
    ## [65] rlang_0.4.10                RCurl_1.98-1.2             
    ## [67] HGNChelper_0.8.1            grid_4.0.2                 
    ## [69] rappdirs_0.3.3              glmGamPoi_1.3.6            
    ## [71] bitops_1.0-6                rmarkdown_2.6              
    ## [73] gtable_0.3.0                codetools_0.2-16           
    ## [75] DBI_1.1.1                   curl_4.3                   
    ## [77] reshape2_1.4.4              R6_2.5.0                   
    ## [79] gridExtra_2.3               knitr_1.31                 
    ## [81] dplyr_1.0.4                 future.apply_1.7.0         
    ## [83] fastmap_1.1.0               bit_4.0.4                  
    ## [85] stringi_1.5.3               parallel_4.0.2             
    ## [87] Rcpp_1.0.6                  sctransform_0.3.2          
    ## [89] vctrs_0.3.6                 dbplyr_2.1.0               
    ## [91] tidyselect_1.1.0            xfun_0.21

</details>

<div id="refs" class="references">

<div id="ref-skene_2016">

1\. Skene, N. & Grant, S. Identification of vulnerable cell types in
major brain disorders using single cell transcriptomes and expression
weighted cell type enrichment. *Frontiers in Neuroscience* (2016).
doi:[10.3389/fnins.2016.00016](https://doi.org/10.3389/fnins.2016.00016)

</div>

<div id="ref-zeisel2015cell">

2\. Zeisel, A. *et al.* Cell types in the mouse cortex and hippocampus
revealed by single-cell rna-seq. *Science* **347,** 1138–1142 (2015).

</div>

</div>
