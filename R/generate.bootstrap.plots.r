#' Generate bootstrap plots
#'
#' \code{generate.bootstrap.plots} takes a genelist and a single cell type transcriptome dataset
#' and generates plots which show how the expression of the genes in the list compares to those
#' in randomly generated gene lists
#'
#' @param sct_data List generated using \code{\link{generate.celltype.data}}
#' @param hits Array of MGI/HGNC gene symbols containing the target gene list.
#' @param bg Array of MGI/HGNC gene symbols containing the background gene list.
#' @param genelistSpecies Either 'mouse' or 'human' depending on whether MGI or HGNC symbols are used for gene lists
#' @param sctSpecies  Either 'mouse' or 'human' depending on whether MGI or HGNC symbols are used for the single cell dataset
#' @param reps Number of random gene lists to generate (default=100 but should be over 10000 for publication quality results)
#' @param annotLevel an integer indicating which level of the annotation to analyse. Default = 1.
#' @param full_results The full output of \code{\link{bootstrap.enrichment.test}} for the same genelist
#' @param listFileName String used as the root for files saved using this function
#' @param savePath Directory where the BootstrapPlots folder should be saved
#' @return Saves a set of pdf files containing graphs. These will be saved with the filename adjusted using the
#' value of listFileName. The files are saved into the 'BootstrapPlot' folder. The files start with one of the following:
#' \itemize{
#'   \item \code{qqplot_noText}: sorts the gene list according to how enriched it is in the relevant celltype. Plots the value in the target list against the mean value in the bootstrapped lists.
#'   \item \code{qqplot_wtGSym}: as above but labels the gene symbols for the highest expressed genes.
#'   \item \code{bootDists}: rather than just showing the mean of the bootstrapped lists, a boxplot shows the distribution of values
#'   \item \code{bootDists_LOG}: shows the bootstrapped distributions with the y-axis shown on a log scale
#' }
#'
#'
#' @examples
#' # Load the single cell data
#' data("ctd")
#'
#' # Set the parameters for the analysis
#' reps=100 # <- Use 100 bootstrap lists so it runs quickly, for publishable analysis use >10000
#'
#' # Load the gene list and get human orthologs
#' data("example_genelist")
#' data("mouse_to_human_homologs")
#' m2h = unique(mouse_to_human_homologs[,c("HGNC.symbol","MGI.symbol")])
#' mouse.hits = unique(m2h[m2h$HGNC.symbol %in% example_genelist,"MGI.symbol"])
#' human.hits = unique(m2h[m2h$HGNC.symbol %in% example_genelist,"HGNC.symbol"])
#' human.bg = unique(m2h$HGNC.symbol)
#' mouse.bg  = unique(m2h$MGI.symbol)
#'
#' # Bootstrap significance testing, without controlling for transcript length and GC content
#' full_results = bootstrap.enrichment.test(sct_data=ctd,hits=mouse.hits,
#'   bg=mouse.bg,reps=reps,annotLevel=1,sctSpecies="mouse",genelistSpecies="mouse")
#'
#' generate.bootstrap.plots(sct_data=ctd,hits=mouse.hits, bg=mouse.bg,
#'      reps=reps,full_results=full_results,listFileName="Example",
#'      genelistSpecies="mouse",sctSpecies="mouse", annotLevel=1)
#' @export
#' @import ggplot2
#' @importFrom reshape2 melt
# @import plyr
generate.bootstrap.plots <- function(sct_data, hits, bg, genelistSpecies="mouse",sctSpecies="mouse", reps, annotLevel=1, full_results=NA, listFileName="",savePath="~/"){
    # Check the arguments
    if(length(full_results)!=3){stop("ERROR: full_results is not valid output from the bootstrap.enrichment.test function")}
    if(sum(!as.character(unique(full_results$results$CellType)) %in% colnames(sct_data[[annotLevel]]$specificity))==length(as.character(unique(full_results$results$CellType)))){
        stop("ERROR: No celltypes in full_results are found in sct_data. Perhaps the wrong annotLevel was used?")
    }
    if(sum(!as.character(unique(full_results$results$CellType)) %in% colnames(sct_data[[annotLevel]]$specificity))>0){
        stop("ERROR: Not all celltypes in full_results are found in sct_data. Perhaps the wrong annotLevel was used?")
    }

    # Add annotLevel to file name tag
    listFileName = sprintf("%s_level%s",listFileName,annotLevel)

    results = full_results$results

    # Check gene lists
    checkedLists = EWCE::check.ewce.genelist.inputs(sct_data,hits,bg,genelistSpecies,sctSpecies)
    hits=checkedLists$hits
    bg=checkedLists$bg
    combinedGenes=unique(c(hits,bg))

    # Get expression data of bootstrapped genes
    signif_res = rownames(results)[results$p<0.05]
    nReps = 1000
    exp_mats = list()
    for(cc in signif_res){
        exp_mats[[cc]] = matrix(0,nrow=nReps,ncol=length(hits))
        rownames(exp_mats[[cc]]) = sprintf("Rep%s",1:nReps)
    }
    for(s in 1:nReps){
        bootstrap_set = sample(combinedGenes,length(hits))
        ValidGenes = rownames(sct_data[[annotLevel]]$specificity)[rownames(sct_data[[annotLevel]]$specificity) %in% bootstrap_set]

        expD = sct_data[[annotLevel]]$specificity[ValidGenes,]
        for(cc in signif_res){
            exp_mats[[cc]][s,] = sort(expD[,cc])
        }
    }


    # Get expression levels of the hit genes
    hit.exp = sct_data[[annotLevel]]$specificity[hits,]

    graph_theme = theme_bw(base_size = 12, base_family = "Helvetica") +
        theme(panel.grid.major = element_line(size = .5, color = "grey"),
              axis.line = element_line(size=.7, color = "black"),legend.position = c(0.75, 0.7), text = element_text(size=14),
              axis.title.x = element_text(vjust = -0.35), axis.title.y = element_text(vjust = 0.6)) + theme(legend.title=element_blank())

    if (!file.exists(sprintf("%s/BootstrapPlots",savePath))){
        dir.create(file.path(savePath, "BootstrapPlots"))
    }

    # Plot the QQ plots
    for(cc in signif_res){
        mean_boot_exp = apply(exp_mats[[cc]],2,mean)
        hit_exp = sort(hit.exp[,cc])
        hit_exp_names = rownames(hit.exp)[order(hit.exp[,cc])]#names(hit_exp)
        dat = data.frame(boot=mean_boot_exp,hit=hit_exp,Gnames=hit_exp_names)
        dat$hit = dat$hit*100
        dat$boot = dat$boot*100
        maxHit = max(dat$hit)
        maxX = max(dat$boot)+0.1*max(dat$boot)
        #if(sub==TRUE){dat$hit=dat$hit/10;dat$boot=dat$boot/10}
        # Plot several variants of the graph
        #basic_graph = ggplot(dat,aes(x=boot,y=hit))+geom_point(size=1)+xlab("Mean Bootstrap Expression")+ylab("Expression in cell type (%)\n") + graph_theme +
        #    geom_abline(intercept = 0, slope = 1, colour = "red")
        basic_graph = ggplot(dat,aes_string(x="boot",y="hit"))+geom_point(size=1)+xlab("Mean Bootstrap Expression")+ylab("Expression in cell type (%)\n") + graph_theme +
            geom_abline(intercept = 0, slope = 1, colour = "red")


        # Plot without text
        pdf(sprintf("%s/BootstrapPlots/qqplot_noText____%s____%s.pdf",savePath,listFileName,cc),width=3.5,height=3.5)
        print(basic_graph)
        dev.off()

        dat$symLab = ifelse(dat$hit>25,sprintf("  %s", dat$Gnames),'')

        #basic_graph = ggplot(dat,aes(x=boot,y=hit))+geom_point(size=2)+xlab("Mean Bootstrap Expression")+ylab("Expression in cell type (%)\n") + graph_theme +
        basic_graph = ggplot(dat,aes_string(x="boot",y="hit"))+geom_point(size=2)+xlab("Mean Bootstrap Expression")+ylab("Expression in cell type (%)\n") + graph_theme +
            geom_abline(intercept = 0, slope = 1, colour = "red")

        # Plot with gene names
        pdf(sprintf("%s/BootstrapPlots/qqplot_wtGSym____%s____%s.pdf",savePath,listFileName,cc),width=3.5,height=3.5)
        print(basic_graph +
                  geom_text(aes_string(label="symLab"),hjust=0,vjust=0,size=3) + xlim(c(0,maxX))
        )
        dev.off()

        # Plot with bootstrap distribution
        melt_boot = melt(exp_mats[[cc]])
        colnames(melt_boot) = c("Rep","Pos","Exp")
        actVals = data.frame(pos=as.factor(1:length(hit_exp)),vals=hit_exp)
        #if(sub==TRUE){melt_boot$Exp=melt_boot$Exp/10; actVals$vals=actVals$vals/10}
        pdf(sprintf("%s/BootstrapPlots/bootDists____%s____%s.pdf",savePath,listFileName,cc),width=3.5,height=3.5)
        melt_boot$Pos = as.factor(melt_boot$Pos)
        #print(ggplot(melt_boot)+geom_boxplot(aes(x=as.factor(Pos),y=Exp),outlier.size=0)+
        #         geom_point(aes(x=pos,y=vals),col="red",data=actVals)+
        print(ggplot(melt_boot)+geom_boxplot(aes_string(x="Pos",y="Exp"),outlier.size=0)+
                  geom_point(aes_string(x="pos",y="vals"),col="red",data=actVals)+
                  ylab("Expression in cell type (%)\n")+
                  xlab("Least specific --> Most specific") + scale_x_discrete(breaks=NULL)+graph_theme)
        dev.off()

        # Plot with LOG bootstrap distribution
        # - First get the ordered gene names
        rownames(dat)=dat$Gnames
        datOrdered = data.frame(GSym=rownames(dat),Pos=1:dim(dat)[1])

        # - Arrange the data frame for plotting
        melt_boot = melt(exp_mats[[cc]])
        colnames(melt_boot) = c("Rep","Pos","Exp")
        melt_boot$Exp = melt_boot$Exp*100
        melt_boot = merge(melt_boot,datOrdered,by="Pos")
        melt_boot$GSym = factor(as.character(melt_boot$GSym),levels=as.character(datOrdered$GSym))

        # - Prepare the values of the list genes to be plotted as red dots
        actVals = data.frame(Pos=as.factor(1:length(hit_exp)),vals=hit_exp*100)
        actVals = merge(actVals,datOrdered,by="Pos")
        actVals$GSym = factor(as.character(actVals$GSym),levels=as.character(datOrdered$GSym))

        # - Determine whether changes are significant
        p = rep(1,max(melt_boot$Pos))
        for(i in 1:max(melt_boot$Pos)){
            p[i] = sum(actVals[actVals$Pos==i,"vals"]<melt_boot[melt_boot$Pos==i,"Exp"])/length(melt_boot[melt_boot$Pos==i,"Exp"])
        }
        ast = rep("*",max(melt_boot$Pos))
        ast[p>0.05] = ""
        actVals = cbind(actVals[order(actVals$Pos),],ast)
        # - Plot the graph!
       # if(sub==TRUE){melt_boot$Exp=melt_boot$Exp/10; actVals$vals=actVals$vals/10}
        wd = 1+length(unique(melt_boot[,4]))*0.175
        pdf(sprintf("%s/BootstrapPlots/bootDists_LOG____%s____%s.pdf",savePath,listFileName,cc),width=wd,height=4)
        #melt_boot$Exp=melt_boot$Exp+0.00000001
        melt_boot = melt_boot[melt_boot$Exp!=0,]
        print(ggplot(melt_boot)+geom_boxplot(aes_string(x="GSym",y="Exp"),outlier.size=0)+graph_theme+
                  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
                  geom_point(aes_string(x="GSym",y="vals"),col="red",data=actVals)+
                  geom_text(aes_string(x="GSym",y="vals",label="ast"),colour="black",data=actVals)+
                  ylab("Expression in cell type (%)\n")+
                  xlab("Least specific --> Most specific")+scale_y_log10()
        )
        dev.off()
    }
}
