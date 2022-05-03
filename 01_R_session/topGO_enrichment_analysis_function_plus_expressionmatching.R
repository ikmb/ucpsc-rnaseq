library(DESeq2)
library(topGO)
library(tidyverse)
library(data.table)
topGO_enrichment <- function(genelist = NULL, 
                             weights = NULL,
                             weight_threshold = 0,
                             statistical_test="Fisher", 
                             algorithm_topGO = "classic",
                             DESeq_result = NULL,
                             match_by_expression = FALSE,
                             gene_background = NULL,
                             ontology_type = "BP",
                             draw_plot = FALSE,
                             output_dir = "output/", 
                             debug_mode = FALSE) {
  # this function is purposed to take a list of genes, which is then analysed by the TOPGO algorithm using GO libraries
  # desired output is a tabular enrichment result
  # also, if the user wishes, a pdf containing a graph of the enrichment hierarchy is written
  # genelist  character vector coontaining gene names. If gene_background is given, genes not in genelist will be ignored.
  # If DESeq_result is given, genes not occuring in the DESeq-assay table will be ignored.
  # weigths are assumed to increasing with importance: Gene A: 5,  Gene B: 2.3: Gene A is more important and will be ranked higher
  # by the Kolmogorov-Smirnov test.
  
  if (debug_mode == TRUE) {
    #this is to intialize for a test run
    genelist <- c("S100A12","ETS1","ACTB")
    weights <- c(2,1,0.8)
    gene_background <- c("VSTM1","TARM1","OSCAR","NDUFA3","TFPT","PRPF31","AC012314.8","CNOT3","LENG1","TMC4","MBOAT7","TSEN34","S100A12","ETS1","ACTB")
  }
  # checks
  stopifnot( "no genelist provided" = !is.null(genelist) )
  stopifnot( "weights and genelist length differ" = is.null(weights) | length(weights) == length(genelist) )
  if (is.null(DESeq_result) & is.null(gene_background)) {
    stop("TopGO needs a gene universe (background). Provide either DESeq_result or gene_background.")
  }
  if ( !is.null(DESeq_result) & !is.null(gene_background) ) {
    message("DESeq_result and gene_background are both input. Using gene_background...")
  }
  library(DESeq2)
  library(topGO)
  library(genefilter)
  library(geneplotter)
  
  # to generate the gene background, prefer gene_background if provided
  # else if no expression matching is wished, use any expressed gene from DESeq2 result
  # else, calculate a background based on similarly expressed genes
  if ( !is.null(gene_background)) {
    gene_background <- gene_background
    backG <- NULL
  } else {
    
    if (match_by_expression == FALSE) {
      DESeq_result <- subset(DESeq_result, DESeq_result$baseMean > 0) 
      gene_background <- DESeq_result$gene
      backG <- NULL
    } else {
      DESeq_result <- subset(DESeq_result, DESeq_result$baseMean > 0) 
      DESeq_result <- subset(DESeq_result, !is.na(DESeq_result$gene) & !DESeq_result$gene == "" )
      rownames(DESeq_result) <- DESeq_result$gene
      gene_background <- rownames(DESeq_result)
      overallBaseMean <- as.matrix(DESeq_result[, "baseMean", drop = F])
      rownames(overallBaseMean) <- rownames(DESeq_result)
      sig_idx <- match(genelist, rownames(overallBaseMean))
      
      backG <- c()
      
      for(i in sig_idx){
        ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
        backG <- c(backG, ind)
        
      }
      
      backG <- unique(backG)
      backG <- rownames(overallBaseMean)[backG]
      backG <- setdiff(backG, genelist)
      
      multidensity( list( 
        all= log2(overallBaseMean[,"baseMean"]) ,
        foreground =log2(overallBaseMean[genelist, "baseMean"]), 
        background =log2(overallBaseMean[backG, "baseMean"])), 
        xlab="log2 mean normalized counts", main = "Matching for enrichment analysis")

    }
    
  }
  
  
  if (!is.null(backG)) {
    inUniverse <- gene_background %in% c(genelist, backG)
    inSelection <- gene_background %in% genelist 
  } else {
    inUniverse <- gene_background %in% gene_background
    inSelection <- gene_background %in% genelist
  }
   # create a named vector with either weights or factor values
  
  if ( !is.null(weights) ) {
    alg <- as.numeric( inSelection[inUniverse] )
    names(alg) <- gene_background[inUniverse]
    alg[genelist] <- weights
  } else {
    alg <- factor( as.integer( inSelection[inUniverse] ) )
    names(alg) <- gene_background[inUniverse]
  }
  
  # function to pick genes which are used when weights are given. 
  
  geneSelFunction <- function(allScore) {
    return(allScore > weight_threshold)
  }
  ## prepare topGO object
  tgd <- new( "topGOdata", ontology = ontology_type, allGenes = alg, 
              geneSel = ifelse( !is.null(weights), geneSelFunction, NULL),
              nodeSize=5,
              annot=annFUN.org, mapping="org.Hs.eg.db", ID = "symbol" )
  
  ## run tests
  #resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
  resultTopGO <- runTest(tgd, algorithm = algorithm_topGO, statistic = statistical_test , scoreOrder = "decreasing") #scoreOrder is handdown for the ks test
  
  ## look at results
  if(length(nodes(graph(tgd))) < 200){
    # tab <- GenTable( tgd, Fisher.elim = resultTopGO.elim, 
    #                  Fisher.classic = resultTopGO.classic,
    #                  orderBy = "Fisher.elim" , topNodes = length(nodes(graph(tgd))))
    tab <- GenTable( tgd, 
                     test.pvalue = resultTopGO,
                     orderBy = "Fisher.elim" , topNodes = length(nodes(graph(tgd))))
  }else{
    # tab <- GenTable( tgd, Fisher.elim = resultTopGO.elim, 
    #                  Fisher.classic = resultTopGO.classic,
    #                  orderBy = "Fisher.elim" , topNodes = 200)
    tab <- GenTable( tgd, 
                     test.pvalue = resultTopGO,
                     orderBy = "Fisher.elim" , topNodes = 200)
    }
  
  
  #printGraph(tgd, resultTopGO.elim, firstSigNodes = 5, fn.prefix = output_dir, useInfo = "all", pdfSW = TRUE)
  if(draw_plot == TRUE) 
    { printGraph(tgd, resultTopGO.classic, firstSigNodes = 15, fn.prefix = output_dir, useInfo = "all", pdfSW = TRUE) }
  
  
  return(tab)
}
topGO_enrichment(debug_mode = TRUE, statistical_test = "ks.ties")
genes.example <- data.table::fread(file = "~/Desktop/example.genes", header = F)$V1
gene.background <- data.table::fread(file = "/home/florian/00_projects/UC_SCZ_RNASeq_paper/uc-rnaseq/00_RawData/genelist.example", header = F)$V1
topGO_enrichment(genelist = genes.example, statistical_test = "fisher", gene_background = gene.background)
res_sorted <- read.table("/home/florian/01_lehre/bioinformatics_module/handson_tutorials/burnhot/output/DESeq2result_genenames_Control_vs_COVID_1.txt", header = TRUE, sep = '\t')
res_sorted <- as.data.table(res_sorted)
res_sorted <- unique(res_sorted, by = "gene")
genelist <- res_sorted[1:100,gene] 
genelist <- genelist[genelist!=""]
weights_res <- res_sorted[res_sorted$gene %in% genelist, stat] %>% abs()
topGO_enrichment(genelist = genelist, DESeq_result = res_sorted, weights =  weights_res, statistical_test = "ks")
topGO_enrichment(genelist = genelist, DESeq_result = res_sorted, weights =  weights_res, statistical_test = "fisher", match_by_expression = TRUE)
