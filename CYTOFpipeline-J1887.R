#Sarah Shin
#Last updated 7/11/23 Draft Run
#R version 4.3.0 (2023-04-21 ucrt)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 19042)


rm(list = ls())
####READ and CLUSTER FUNCTIONS####
returnfcs <- function(FDR_cutoff=.05,
                      metaDataFile='~/',
                      panelDataFile='~/',
                      dataDirectory='~/',
                      shape_conditions=NULL,
                      color_conditions=NULL){
  #This function generates an fcs file, subtype_markers, colors and shapes for clustering 
  require(scales);require(readxl);require(dplyr);require(flowCore)
  ##directory and metadatafile checking
  if(!dir.exists(dataDirectory)) {stop('ERR: cannot find data directory')}
  if(!file.exists(metaDataFile)) {stop('ERR: cannot find metadata.xlsx or .csv file')}
  ##readin metadata and clean
  ifelse(grepl(metaDataFile,pattern='.xls'),md <- read_excel(metaDataFile),md <- read.csv(metaDataFile,header = TRUE))#must be in xl format or csv
  md$response <- factor(md$response);md$survival<-factor(md$survival);md$timepoint <- factor(md$timepoint)
  rownames(md) = md$sample_id;md$sample_id <- md$sample_id
  #Make sure all files in metadata present in datadirectory
  if(!all(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])){
    print(paste('ERR: not all filenames in metadata present in data folder - missing',md$file_name[!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])],'Subsetting...'))
    md <- md[-c(!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])),]
  }
  ##Define shapes for conditions
  if(is.null(shape_conditions)){shape_conditions <- c(0:25)[1:length(levels(md$response))]}#can specify as long as number is same
  if(length(shape_conditions)!=length(levels(md$response))){stop(paste0('ERR no. shapes specified is less than no. of conditions (',length(levels(md$response)),')'))}
  names(shape_conditions) <- levels(md$response)
  ## Define colors for the conditions
  if(is.null(color_conditions)){color_conditions <- hue_pal()(length(levels(md$response)))}#can specify as long as number is same
  if(length(color_conditions)!=length(levels(md$response))){stop(paste0('ERR no. shapes specified is less than no. of conditions (',length(levels(md$response)),')'))}
  ## read fcs
  fcs_raw <- read.flowSet(md$file_name, path = dataDirectory, transformation = FALSE, truncate_max_range = FALSE)
  sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))
  panel <- read_excel(panelDataFile)
  head(data.frame(panel))
  ## Replace problematic characters
  panel$Metal <- gsub('-', '_', panel$Metal)
  panel_fcs <- pData(parameters(fcs_raw[[1]]))
  panel_fcs$desc <- gsub('-', '_', panel_fcs$desc)
  panel_fcs$desc[is.na(panel_fcs$desc)] <- paste0('NA_',which(is.na(panel_fcs$desc))) #was labelled 'todo'(mapping based on isotope for now just getting rid of NA keeping rownum) 
  # use panel$Antigen to fix description in panel_fcs
  # use metal+isotope as mapping between panel from xlsx and panel from the fcs files
  rownames(panel_fcs) = panel_fcs$name
  panel_fcs[paste0(panel$Metal,panel$Isotope,'Di'),2] <- panel$Antigen
  ## Replace paramater data in flowSet
  pData(parameters(fcs_raw[[1]])) <- panel_fcs
  ## Define variables indicating marker types
  subtype_markers <- panel$Antigen[panel$Subtype == 1]
  functional_markers <- panel$Antigen[panel$Functional == 1]
  if(!all(subtype_markers %in% panel_fcs$desc)){stop('ERR: Not all subtype_markers in panel_fcs$desc (isotopes)')}
  if(!all(functional_markers %in% panel_fcs$desc)){stop('ERR: Not all functional_markers in panel_fcs$desc (isotopes)')}
  ## arcsinh transformation and column subsetting
  fcs <- fsApply(fcs_raw, function(x, cofactor = 5){
    colnames(x) <- panel_fcs$desc
    expr <- exprs(x)
    expr <- asinh(expr[, union(subtype_markers,functional_markers)] / cofactor)
    exprs(x) <- expr
    x
  })
  return(list('fcs'=fcs,
              'subtype_markers'=subtype_markers,
              'functional_markers'=functional_markers,
              'shape_conditions'=shape_conditions,
              'color_conditions'=color_conditions,
              'sample_ids'=sample_ids,
              'meta_data'=md))
}

clusterfcs <- function(fcs=output$fcs,
                       subtype_markers = output$subtype_markers,
                       seed=1234,plottitle='consensus_plots',
                       numclusters=40){
  ## Cell population identification with FlowSOM and ConsensusClusterPlus
  require(dplyr);require(FlowSOM);require(ConsensusClusterPlus)
  set.seed(seed)
  som <- ReadInput(fcs, transform = FALSE, scale = FALSE) %>% BuildSOM(colsToUse = subtype_markers)
  ## Get the cell clustering into 100 SOM codes
  cell_clustering_som <- som$map$mapping[,1]
  ## Metaclustering into numclusters with ConsensusClusterPlus
  codes <- som$map$codes
  mc <- ConsensusClusterPlus(t(codes), maxK = numclusters, reps = 100,
                             pItem = 0.9, pFeature = 1, title = plottitle, 
                             plot = "png", clusterAlg = "hc", 
                             innerLinkage = "average", finalLinkage = "average",
                             distance = "euclidean", seed = 1234)
  
  ## Get cluster ids for each cell
  code_clustering <- mc[[numclusters]]$consensusClass#metaclusters consensus
  cell_clustering <- code_clustering[cell_clustering_som]#cell clustering from som
  return(list('code_clustering'=code_clustering,'cell_clustering'=cell_clustering,'metaclusters'=mc))
}



####DIAGNOSTIC HEATMAP FUNCTIONS ####
plot_clustering_heatmap_wrapper <- function(fcs, cell_clustering, nclusters=40,
                                            cluster_merging = NULL, 
                                            subtype_markers,
                                            clusterMergeFile=NULL,
                                            fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales)
  ## Will output the heatmap object and print it 
  color_clusters <-c(stepped2(20),stepped3(20),rev(cubicl(20)))
  #get expression
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.01, 0.99))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,subtype_markers]
  ## Calculate the mean expression##################################################
  pdf(fileName, width=8, height=11) 
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Colors for the heatmap
  
  color_heat <- colorRampPalette(brewer.pal(n = 9, name = "YlOrBr"))(100)
  #legend_breaks = seq(from = 0, to = 1, by = 0.2)
  labels_row <- paste0(expr01_mean$cell_clustering, " (", clustering_prop ,
                       "%)")
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  ## Annotation for the merged clusters
  
  if(!is.null(clusterMergeFile)){
    ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging <- read_excel(clusterMergeFile),cluster_merging <- read.csv(clusterMergeFile,header = TRUE))
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$Merged <- cluster_merging$new_cluster
    color_clusters2 <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters2) <- levels(cluster_merging$new_cluster)
    annotation_colors$Merged <- color_clusters2
  }
  
  p <- pheatmap(expr_heat, color = rev(colorRampPalette(brewer.pal(n = 11, name = "RdYlBu"))(100)), 
                cluster_cols = T,
                cluster_rows = T, 
                labels_row = labels_row,
                #scale="column",
                display_numbers = FALSE, number_color = "black",
                fontsize = 9, fontsize_number = 6,  
                #legend_breaks = legend_breaks,
                annotation_row = annotation_row, 
                annotation_colors = annotation_colors,
                cellwidth = 8,
                cellheight = 8
  )
  dev.off() 
  print('Colors:')
  print(color_clusters)
  print(p);return(p)
}


plot_clustering_heatmap_wrapper2 <- function(fcs, cell_clustering, nclusters=40,
                                             color_clusters=clustercolors,
                                             subtype_markers,
                                             fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales)
  ## Will output the heatmap object and print it 
  color_clusters=color_clusters
  #get expression
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.01, 0.99))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,subtype_markers]
  ## Calculate the mean expression##################################################
  pdf(fileName, width=8, height=11) 
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Colors for the heatmap
  
  #legend_breaks = seq(from = 0, to = 1, by = 0.2)
  #labels_row <- expr01_mean$cell_clustering
  
  labels_row <- paste0(expr01_mean$cell_clustering, " ")
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  p <- pheatmap(expr_heat, 
                color = c(rep(magma(100)[1],25),magma(100)[1:100]), 
                cluster_cols = T,
                cluster_rows = F, 
                labels_row = labels_row,
                #scale="column",
                display_numbers = F, 
                number_color = "black",
                fontsize = 9, fontsize_number = 6,  
                #legend_breaks = legend_breaks,
                annotation_row = annotation_row, 
                annotation_colors = annotation_colors,
                cellwidth = 8,
                cellheight = 8,
                border_color = "black",
                annotation_legend = F
  )
  dev.off() 
  print('Colors:')
  print(color_clusters)
  print(p);return(p)
}


####DIAGNOSTICS####
makeDiagnosticPlots = function(exprData, 
                               md = output$meta_data,
                               sample_ids = output$sample_ids,
                               fcs = output$fcs,
                               samplevels = samplevels,
                               subtype_markers = output$subtype_markers,
                               color_conditions = clustercolors,
                               shape_conditions = c(1:13),
                               fileName = 'diagnostics.pdf', 
                               tit = '', 
                               fun = mean)
{
  pdf(file = fileName)
  
  ## Spot check - number of cells per sample
  cell_table <- table(sample_ids)
  ggdf <- data.frame(sample_id = names(cell_table), 
                     cell_counts = as.numeric(cell_table))
  ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels)
  mm <- match(ggdf$sample_id, md$sample_id)
  ggdf$condition <- md$condition[mm]
  print(ggplot(ggdf, aes(x = sample_id, y = cell_counts, fill = condition)) + 
          geom_bar(stat = 'identity') + 
          geom_text(aes(label = cell_counts), hjust = 0.5, vjust = -0.5, size = 2.5) + 
          theme_bw() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +  
          scale_fill_manual(values = color_conditions, drop = FALSE) + 
          scale_x_discrete(drop = FALSE))
  
  dev.off()
}


####CLUSTER HISTO####

plot_clustering_distr_wrapper <- function(expr = expr, 
                                          cell_clustering){
  # Calculate the median expression
  cell_clustering <- factor(cell_clustering)
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  # Sort the cell clusters with hierarchical clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  # Calculate cluster frequencies
  freq_clust <- table(cell_clustering)
  freq_clust <- round(as.numeric(freq_clust)/sum(freq_clust)*100, 2)
  cell_clustering <- factor(cell_clustering,
                            labels = levels(cell_clustering))
  ### Data organized per cluster
  ggd <- melt(data.frame(cluster = cell_clustering, expr),
              id.vars = "cluster", value.name = "expression",
              variable.name = "antigen")
  ggd$antigen <- factor(ggd$antigen, levels = colnames(expr))
  ggd$reference <- "no"
  ### The reference data
  ggd_bg <- ggd
  ggd_bg$cluster <- "reference"
  ggd_bg$reference <- "yes"
  
  ggd_plot <- rbind(ggd, ggd_bg)
  ggd_plot$cluster <- factor(ggd_plot$cluster,
                             levels = c(levels(cell_clustering)[rev(cluster_rows$order)], "reference"))
  
  ggplot() +
    geom_density_ridges(data = ggd_plot, aes(x = expression, y = cluster,
                                             color = reference, fill = reference), alpha = 0.3) +
    facet_wrap( ~ antigen, scales = "free_x", nrow = 2) +
    theme_ridges() +
    theme(axis.text = element_text(size = 7),  
          strip.text = element_text(size = 7), legend.position = "none")
  
}

####DIAGNOSTICS####
makeDiagnosticPlots = function(exprData, 
                               md = output$meta_data,
                               sample_ids = output$sample_ids,
                               fcs = output$fcs,
                               samplevels = samplevels,
                               subtype_markers = output$subtype_markers,
                               color_conditions = clustercolors,
                               shape_conditions = c(1:13),
                               fileName = 'diagnostics.pdf', 
                               tit = '', 
                               fun = mean)
{
  pdf(file = fileName)
  
  ## Spot check - number of cells per sample
  cell_table <- table(sample_ids)
  ggdf <- data.frame(sample_id = names(cell_table), 
                     cell_counts = as.numeric(cell_table))
  ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels)
  mm <- match(ggdf$sample_id, md$sample_id)
  ggdf$condition <- md$condition[mm]
  print(ggplot(ggdf, aes(x = sample_id, y = cell_counts, fill = condition)) + 
          geom_bar(stat = 'identity') + 
          geom_text(aes(label = cell_counts), hjust = 0.5, vjust = -0.5, size = 2.5) + 
          theme_bw() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +  
          scale_fill_manual(values = color_conditions, drop = FALSE) + 
          scale_x_discrete(drop = FALSE))
  
  dev.off()
}

####UMAP + HEXBIN####
do_umap <- function(fcs,subtype_markers,sample_ids,cell_clustering,metadata,
                    clusterMergeFile='~/Desktop/ViralHCC/ViralHCC_merging.xlsx',
                    seed = 1234, ncells=2000,sample_subset=NULL, neighbors=10){
  require(umap);require(flowCore);require(readxl)
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Create vector to later find and skip duplicates
  dups <- duplicated(expr[, subtype_markers])
  dups <- which(!(dups))## Find and skip duplicates
  ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging <- read_excel(clusterMergeFile),cluster_merging <- read.csv(clusterMergeFile,header = TRUE))
  ## New clustering1m
  mm <- match(cell_clustering, cluster_merging$original_cluster)
  cell_clustering1m <- cluster_merging$new_cluster[mm]
  ## Create a data frame of sample_ids and cell_clustering1m
  dtf<-data.frame(ids=sample_ids,type=cell_clustering1m)
  #dtf$B<- dtf$type!="B"#add a column that indicates whether the cell is a B cell or not; TRUE is non-B
  ##Why exclude B CELLS?
  ## WE HAVE NO B CELLS bc dtf$B depends on dtf$type depends on cellclustering1m which is just numbers 1:40 so..?
  #should it be the named parts cluster in merge file corresponding to it like if 30 -> grepl(cluster_merging[30,3],pattern='^B')??
  #Im blocking this out till we know why we have to do this
  ## Create subset columns without B cells (ONLY to generate the correct # of columns in inds2 object)
  #sample_ids2 <- dtf$ids[dtf$type!="B"] #sampleids without B cells
  #cell_clustering1m2 <- dtf$type[dtf$type!="B"] #clusters without B cells
  ## Data subsampling: create indices by sample
  inds <- split(1:length(sample_ids), sample_ids) #to get original indexes belonging to each cluster
  #inds2 <- split(1:length(sample_ids2), sample_ids2) #to fill in the original indexes that do not have B cells
  samplenames <- names(inds) #create a name vector of the files
  #FOR THIS TO WORK MUST BE IN FORMAT PRE/POST Tx
  # for (i in 1:(length(samplenames)/2)){#1:15 was here because ids in dtf was 30 factors and could be divided into 0 and 6wk for each so changed it
  #   templength <- length(inds2[[i]])
  #   inds2[[i]] <- inds[[i]][dtf$B[dtf$ids==samplenames[i]]] #subsets the "B cell or not a B cell" column for each sample by the name
  #   inds2[[i]] <- inds2[[i]][1:templength]
  # }
  custom.settings = umap.defaults
  custom.settings$seed = seed
  custom.settings$n.neighbors = neighbors
  ####umapindex generation####
  #umap ncells = table of sample ids with how many to downsample to by default col = id, row = ncells
  #sample ids = chr [1:2851129] "4927_0wk" "4927_0wk" "4927_0wk" "4927_0wk" ...
  #^ from ## Generate sample IDs corresponding to each cell in the 'expr' matrix sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))
  #can subset sample_ids and rerun umap 
  #can do timepoint or patient number or pre/post if you know corresp sample ids
  #sample_subset = '02_0wk' or c('02_0wk','2973_0wk') for example picking all 0 wk ones or use regex sample_ids[(grepl(sample_ids,pattern = '0wk'))]
  ifelse(is.null(sample_subset),
         umap_ncells <- pmin(table(sample_ids), ncells),
         umap_ncells <- pmin(table(sample_ids), ncells)[sample_subset]
  )
  if(!is.null(sample_subset)){inds <- inds[sample_subset]}
  umap_inds <- lapply(names(inds), function(i){
    s <- sample(inds[[i]], umap_ncells[i], replace = FALSE)
    intersect(s, dups)
  })
  set.seed(seed)
  umap_inds <- unlist(umap_inds)
  umap_out <- umap(expr[umap_inds, subtype_markers], config = custom.settings, method = 'naive')
  umapRes2D = data.frame(umap1 = umap_out$layout[, 1], umap2 = umap_out$layout[, 2], 
                         expr[umap_inds, subtype_markers],
                         sample_id = sample_ids[umap_inds], cell_clustering = factor(cell_clustering1m[umap_inds]), check.names = FALSE)
  
  #exclude any unassigned cluster post umap if needed--this has to be done by looking at the two columns 
  #metadata$sampleid is just a number in this metadatafile so to make unique ones combine samp_id col with timepoint (0wk)
  #to get in format of umapres2d$sample_id which looks like "02_0wk" do:
  return(umapRes2D)
}



plotUmap <- function(umapRes,
                     seed=1234,
                     neighbors=10,
                     midpoint,
                     color_clusters,
                     code_clustering,subtype_markers=NULL)
{require(umap);require(ggplot2);require(viridis);require(ggrepel)
  custom.settings = umap.defaults
  custom.settings$seed = seed
  custom.settings$n.neighbors = neighbors
  ggp <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = cell_clustering)) +
    geom_point(size = 0.1) +
    theme_bw() +
    
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key.size=unit(1, "lines")
    ) +
    
    scale_color_manual(values = color_clusters, name="CLUSTERS") +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))
  
  print(ggp)
  #other options
  print(ggp + facet_wrap(~ timepoint)+ggtitle('PRE vs. POST'))
  
  print(ggp + facet_wrap(~ survival)+ggtitle('SHORT vs. LONG'))

  print(ggp + facet_wrap(~ survival ~timepoint))
  
  print(ggp + facet_wrap(~ batch)+ggtitle('BATCH'))
  
  print(ggp + facet_wrap(~patient_id) + ggtitle('PATIENT'))
  
  ggp2 <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = batch)) +
    geom_point(size = 0.1) +
    theme_bw() +
    
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key.size=unit(0.5, "lines")
    ) +
    
    scale_color_manual(values = color_clusters, name="BATCH") +
    guides(color = guide_legend(override.aes = list(size = 2), ncol = 2))
  
  print(ggp2)
  
  #can specify which markers to display
  if(!is.null(subtype_markers)){
    for(i in subtype_markers)
    {
      ggp <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = umapRes[,i])) +
        geom_point(size = 1) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        scale_color_gradient2(i, low="dark blue",mid="white",high="dark red", midpoint = mean(unlist(umapRes[,i])))
      print(ggp)
    }
  }
}





####DS ANALYSIS####
de_wrapper <- function(expr_mean, md, model = "lmer", formula, K, contrast_names){
  ## Define K by subsetting the matrix, e.g. K[1,] or K[2,]
  ## Define contrast_names by selecting, e.g. contrast_names[1] 
  ## Fit LMM or LM for each marker separately
  fit_gaussian <- lapply(1:nrow(expr_mean), function(i){
    data_tmp <- data.frame(y = as.numeric(expr_mean[i, md$sample_id]), md)
    switch(model,
           lmer = {
             fit_tmp <- lmer(formula, data = data_tmp)
           },
           lm = {
             fit_tmp <- lm(formula, data = data_tmp)
           })
    ## Fit contrasts one by one
    contr_tmp<-glht(fit_tmp,linfct = t(matrix(K)))
    summ_tmp<-summary(contr_tmp)
    pval <- summ_tmp$test$pfunction()
    return(pval)
  })
  pvals <- do.call(rbind, fit_gaussian)
  colnames(pvals) <- paste0("pval_", contrast_names)
  rownames(pvals) <- rownames(expr_mean)
  ## Adjust the p-values
  adjp <- apply(pvals, 2, p.adjust, method = "BH")
  colnames(adjp) <- paste0("adjp_", contrast_names)
  return(list(pvals = pvals, adjp = adjp))
}

####ADDITIONAL LIBRARY LOADING####
library(reshape2);library(diffcyt);library(ggpubr);library(readxl)
library(randomcoloR);library(pals);library(ggplot2);library(colorspace);
library(dplyr);library(scales);library(tibble)
library(lme4);library(multcomp);library(mvtnorm); library(flowCore)
library(umap); library(matrixStats); library(premessa); library(limma);
library(ComplexHeatmap);library(circlize); library(stringr); library(edgeR)
library(corrplot);library(ggridges); library(Hmisc); library(pheatmap); 
library(RColorBrewer);

#==========DATA LOADING AND CLUSTERING===========
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
workd<-getwd()

#if it's needed to load up saved output RDS
output<-readRDS('backup_output.rds')

#read and cluster

output <- returnfcs(metaDataFile = paste0(workd,"/Config","/metadata.xlsx"),
                    panelDataFile = paste0(workd,"/Config","/panel.xlsx"),
                    dataDirectory = paste0(workd,'/Data'))

output[8:10] <- clusterfcs(numclusters=30)
names(output)[8:10] <- c('code_clustering','cell_clustering','metaclusters')


#assign clusters

clusterMergeFile = paste0(workd,"/Config",'/merge.xlsx')

cluster_merging <- read_xlsx(clusterMergeFile)

clusterlevels=c(1:30)
clusterlevels = c("B_I",
                  "B_II",
                  "B_III",
                  "DNT",
                  "DNT_Acv_I",
                  "DNT_Acv_II",
                  "Myeloid",
                  "NK",
                  "NK_Acv",
                  "TcCM",
                  "TcEFF_Acv_I",
                  "TcEFF_Acv_II",
                  "TcEM",
                  "TcN",
                  "ThCM",
                  "ThCM_Acv",
                  "ThEM_Acv_I",
                  "ThEM_Acv_II",
                  "ThEM_Acv_III",
                  "ThN",
                  "Treg",
                  "UA")

samplevels <- c("E1_BC1","E1_BC10","E1_BC11","E1_BC12","E1_BC13","E1_BC14","E1_BC2",
                "E1_BC3","E1_BC4","E1_BC5","E1_BC6","E1_BC7","E1_BC8","E1_BC9",
                "E2_BC1","E2_BC10","E2_BC11","E2_BC12","E2_BC13","E2_BC2",
                "E2_BC3","E2_BC4","E2_BC5","E2_BC6","E2_BC7","E2_BC8","E2_BC9",
                "E3_BC1","E3_BC10","E3_BC11","E3_BC12","E3_BC13","E3_BC14",
                "E3_BC2","E3_BC3","E3_BC4","E3_BC5","E3_BC6","E3_BC7","E3_BC8","E3_BC9")

ptlevels <-c("002_C1D1","164_C1D1","164_C3D1","254_C1D1","254_C3D1","258_C1D1","002_C3D1",
             "205_C1D1","205_C3D1","001_C1D1","153_C1D1","153_C3D1","159_C1D1","159_C3D1",
             "102_C1D1","165_C1D1","165_C3D1","255_C1D1","255_C3D1","102_C3D1","206_C1D1","206_C3D1",
             "151_C3D1","154_C1D1","154_C3D1","160_C1D1","160_C3D1","103_C1D1","167_C1D1","251_C1D1",
             "256_C3D1","257_C1D1","257_C3D1","103_C3D1","207_C1D1","207_C3D1","152_C1D1","155_C1D1",
             "158_C1D1","162_C3D1","163_C1D1")

caselevels <- c("002","164","254","258","205","001","153","159",
              "102","165","255","206","151","154","160","103","167",
              "251","256","257","207","152","155","158","162","163")


condlevels=c("WT","Mutant")

timepointlevels=c("C1D1", "C3D1")

timegrouplevels = c("PreTx", "PostTx")

responselevels=c("PD", "SD", "PR")

benefitlevels=c("Non-Benefit", "Benefit")

survivallevels=c("Long", "Short")

batchlevels=c("1","2","3")

clustercolors <- as.character(c(cols25(n=25),alphabet(n=19)))

mm1 <- match(output$cell_clustering, cluster_merging$original_cluster)

cell_clustering1m <- cluster_merging$new_cluster[mm1]

output$cell_clustering1m <- cell_clustering1m

#metacluster heatmap
plot_clustering_heatmap_wrapper(fcs=output$fcs,
                                cell_clustering = output$cell_clustering, 
                                subtype_markers=output$subtype_markers,
                                clusterMergeFile = clusterMergeFile,
                                fileName = 'clusteringheatmap2.pdf');dev.off()

plot_clustering_heatmap_wrapper2(fcs=output$fcs,
                                 color_clusters = clustercolors,
                                 cell_clustering = factor(output$cell_clustering1m, levels=clusterlevels), 
                                 subtype_markers=output$subtype_markers,
                                 fileName = 'clusteringheatmap_merged.pdf');dev.off()

#save output
saveRDS(output, file = "backup_output.rds")

## DIAGNOSTICS
expr <- fsApply(output$fcs, exprs) #create expression matrix
rng <- colQuantiles(expr, probs = c(0.01, 0.99))
expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1])) #scaling 0-1
expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1 # need this expr object later 

dgplots<-makeDiagnosticPlots(expr01, 
                             md = output$meta_data,
                             sample_ids = output$sample_ids,
                             samplevels = samplevels,
                             fcs = output$fcs,
                             subtype_markers = output$subtype_markers,
                             color_conditions = clustercolors,
                             shape_conditions = c(1:13),
                             fileName = 'diagnostics.pdf')
                             dev.off()

## MDS PLOT
expr_mean_sample_tbl <- data.frame(sample_id = output$sample_ids, expr01) %>%
  group_by(sample_id) %>%  summarize_all(funs(mean))

expr_mean_sample <- t(expr_mean_sample_tbl[, -1])
colnames(expr_mean_sample) <- expr_mean_sample_tbl$sample_id

mds <- plotMDS(expr_mean_sample, plot = FALSE)
ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y,
                   sample_id = colnames(expr_mean_sample))
ggdf$sample_id <- factor(ggdf$sample_id, levels = samplevels)
mm <- match(ggdf$sample_id, output$meta_data$sample_id)
ggdf$condition <- output$meta_data$condition[mm]
ggdf$batch <- output$meta_data$batch[mm]
ggdf$timepoint <- output$meta_data$timepoint[mm]
ggdf$timegroup <- output$meta_data$timegroup[mm]
ggdf$response <- output$meta_data$response[mm]
ggdf$benefit <- output$meta_data$benefit[mm]
ggdf$survival <- output$meta_data$survival[mm]
ggdf$case <- output$meta_data$case[mm]
ggdf$patient_id <- output$meta_data$patient_id[mm]


mdsggp <- ggplot(ggdf, aes(x = MDS1, y = MDS2, color = condition)) +
  geom_point(size = 2, alpha = 0.8) +
  #geom_label_repel(aes(label = sample_id)) +
  theme_bw() +
  theme(axis.text = element_text(color="black")) +
  #ylim(-0.1, 0.08)+
  #xlim(-0.12, 0.12)+
  scale_color_manual(values = c("blue","yellow","green")) +
  scale_shape_manual(values = c(2,0,1)) +
  coord_fixed()
pdf("plot_condition_mds.pdf", width=4, height=5);mdsggp;dev.off()




## UMAP
umapRes <- do_umap(fcs=output$fcs,subtype_markers = output$subtype_markers,
                   sample_ids = output$sample_ids,cell_clustering = output$cell_clustering, metadata=output$metadata,
                   clusterMergeFile=clusterMergeFile,
                   seed = 1234, ncells=100,sample_subset=NULL)

mm <- match(as.character(umapRes$sample_id), as.character(output[["meta_data"]]$sample_id))
umapRes$sample_id <- factor(output[["meta_data"]]$sample_id[mm], levels=samplevels)
umapRes$patient_id <- factor(output[["meta_data"]]$patient_id[mm], levels=ptlevels)
umapRes$case <- factor(output[["meta_data"]]$case[mm], levels = caselevels)
umapRes$timepoint <- factor(output[["meta_data"]]$timepoint[mm], levels=timepointlevels)
umapRes$timegroup <- factor(output[["meta_data"]]$timegroup[mm], levels=timegrouplevels)
umapRes$response <- factor(output[["meta_data"]]$response[mm], levels=responselevels)
umapRes$benefit <- factor(output[["meta_data"]]$benefit[mm], levels=benefitlevels)
umapRes$condition <- factor(output[["meta_data"]]$condition[mm], levels=condlevels)
umapRes$survival <- factor(output[["meta_data"]]$survival[mm], levels=survivallevels)
umapRes$batch <- factor(output[["meta_data"]]$batch[mm], levels=batchlevels)


pdf('plot_umaps.pdf',width=10,height=10)
plotUmap(umapRes = umapRes,
         code_clustering=cell_clustering1m,
         color_clusters = clustercolors)

dev.off()

saveRDS(umapRes, file="backup_umap.rds")
umapRes<-readRDS('backup_umap.rds')

#==========DIFFERENTIAL ABUNDANCE PLOTS===========
counts_table <- table(output$cell_clustering1m, output$sample_ids)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)


ggdf <- melt(data.frame(cluster = rownames(props), props, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "sample_id")
ggdf$sample_id <- factor(ggdf$sample_id)
ggdf$cluster <- factor(ggdf$cluster, levels=clusterlevels)
ggdf$timepoint <- factor(output$meta_data$timepoint[match(ggdf$sample_id,output$meta_data$sample_id)], levels = timepointlevels)
ggdf$timegroup <- factor(output$meta_data$timegroup[match(ggdf$sample_id,output$meta_data$sample_id)], levels = timegrouplevels)
ggdf$response <- factor(output$meta_data$response[match(ggdf$sample_id,output$meta_data$sample_id)], levels = responselevels)
ggdf$survival <- factor(output$meta_data$survival[match(ggdf$sample_id,output$meta_data$sample_id)], levels = survivallevels)
ggdf$case <- factor(output$meta_data$case[match(ggdf$sample_id,output$meta_data$sample_id)], levels = caselevels)
ggdf$benefit <- factor(output$meta_data$benefit[match(ggdf$sample_id,output$meta_data$sample_id)], levels = benefitlevels)
ggdf$patient_id <- factor(output$meta_data$patient_id[match(ggdf$sample_id,output$meta_data$sample_id)], levels = ptlevels)
ggdf$condition <- factor(output$meta_data$condition[match(ggdf$sample_id,output$meta_data$sample_id)], levels = condlevels)
ggdf$batch <- factor(output$meta_data$batch[match(ggdf$sample_id,output$meta_data$sample_id)], levels = batchlevels)
ggdf$wbc <- factor(output$meta_data$WBC[match(ggdf$sample_id,output$meta_data$sample_id)])
ggdf$cond_tp <- paste(ggdf$condition, ggdf$timepoint, sep="_")
ggdf$cond_tp <- factor(ggdf$cond_tp, levels = c("WT_C1D1", "WT_C3D1", "Mutant_C1D1", "Mutant_C3D1"))
ggdf$resp_tp <- paste(ggdf$response, ggdf$timepoint, sep="_")
ggdf$resp_tp <- factor(ggdf$resp_tp, levels = c("PD_C1D1", "PD_C3D1", "SD_C1D1", "SD_C3D1", "PR_C1D1", "PR_C3D1"))
ggdf$surv_tp <- paste(ggdf$survival, ggdf$timepoint, sep="_")
ggdf$surv_tp <- factor(ggdf$surv_tp, levels = c("Short_C1D1", "Short_C3D1", "Long_C1D1", "Long_C3D1"))

write.csv(ggdf,"props.csv")


boxplot <- ggplot(ggdf[ggdf$timepoint=="C1D1",], aes(x=survival, y=proportion, fill=survival))+
  geom_boxplot(outlier.size=0, lwd=0.25, color="black")+
  geom_jitter(width=0)+
  facet_wrap(~cluster,ncol=6,scales="free")+
  ylab("% of CD45+ Cells")+
  ggtitle("Pre-Tx by Survival")+
  #scale_fill_manual(values=c("#c93c3c","#3cc95b"))+
  theme(axis.text.x = element_text(size=10, angle=30, vjust=1, hjust=1, color="black"),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.5, color="black"),
        axis.line.y = element_line(size=0.5, color="black"),
        axis.title.y = element_text(size=12, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=7.5, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white")
  )+
  stat_compare_means(ref.group = "Short",
                     paired=F,
                     method="t.test",
                     label = "p.format",
                     label.x = 1.4,
                     label.y.npc = "top",
                     hide.ns = F,
                     size=1.5)
pdf("abundance_box_MERGED_prop_pretx_survival_stats+.pdf", width=11, height=16);boxplot;dev.off()


resp_tp_comparisons <- list(c("PD_C1D1", "PD_C3D1"), c("SD_C1D1", "SD_C3D1"), c("PR_C1D1", "PR_C3D1"))
surv_tp_comparisons <- list(c("Short_C1D1", "Short_C3D1"), c("Long_C1D1", "Long_C3D1"))
cond_tp_comparisons <- list(c("WT_C1D1", "WT_C3D1"), c("Mutant_C1D1", "Mutant_C3D1"))

ggp <- ggplot(ggdf, aes(x=cond_tp, y=proportion, group=case, col=case))+
  facet_wrap(~cluster, scales='free',ncol=6)+
  #scale_color_manual(values=c("#c93c3c","#3cc95b"))+
  scale_shape_manual(values=c(0:40))+
  geom_point()+
  geom_line(show.legend = F)+
  ggtitle("Condition (Proportion)")+
  theme(axis.text.x = element_text(size=8, angle = 45, hjust=1, vjust=1, color="black"),
        axis.title.x = element_blank(),
        axis.ticks=element_line(color="black"),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=10),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(1,'lines'),
        legend.text = element_text(size=10),
        legend.key = element_rect(fill="white"))
pdf("lineplot_prop_all_condition+timepoint.pdf", width=11, height=8);ggp;dev.off()

                       
ggp <- ggplot(ggdf, aes(x=timepoint, y=proportion, group=case, col=condition))+
  facet_wrap(~cluster, scales='free',ncol=6)+
  #scale_color_manual(values=c("#c93c3c","#3cc95b"))+
  scale_shape_manual(values=c(0:24))+
  geom_point()+
  geom_line(show.legend = F)+
  ggtitle("Timepoint by Condition")+
theme(axis.text.x = element_text(size=8, angle = 45, hjust=1, vjust=1, color="black"),
      axis.title.x = element_blank(),
      axis.ticks=element_line(color="black"),
      axis.line.x = element_line(size=0.5),
      axis.line.y = element_line(size=0.5),
      axis.text.y = element_text(color="black"),
      axis.title.y = element_blank(),
      strip.background = element_rect(fill=NA),
      strip.text = element_text(size=10),
      panel.background = element_rect(fill="white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.key.size = unit(1,'lines'),
      legend.text = element_text(size=10),
      legend.key = element_rect(fill="white")+
        stat_compare_means(comparisons = my_comparisons,
                           paired=F,
                           method="t.test",
                           label = "p.format",
                           label.x = 1.4,
                           label.y.npc = "top",
                           hide.ns = F,
                           size=1.5))
pdf("lineplot_prop_all_condition_stats+.pdf", width=11, height=8);ggp;dev.off()

ggp <- ggplot(ggdf, aes(x=timepoint, y=proportion, group=case))+
  facet_wrap(~cluster, scales='free',ncol=6)+
  #scale_color_manual(values=c("#c93c3c","#3cc95b"))+
  scale_shape_manual(values=c(0:24))+
  geom_point(size = 2, aes(color = combresp))+
  geom_line(show.legend = F)+
  ggtitle("xn")+
  scale_x_discrete()+
theme(axis.text.x = element_text(size=8, angle = 45, hjust=1, vjust=1, color="black"),
      axis.title.x = element_blank(),
      axis.ticks=element_line(color="black"),
      axis.line.x = element_line(size=0.5),
      axis.line.y = element_line(size=0.5),
      axis.text.y = element_text(color="black"),
      axis.title.y = element_blank(),
      strip.background = element_rect(fill=NA),
      strip.text = element_text(size=10),
      panel.background = element_rect(fill="white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.key.size = unit(1,'lines'),
      legend.text = element_text(size=10),
      legend.key = element_rect(fill="white"))
pdf("lineplot_prop_combresp_timepoint.pdf", width=11, height=8);ggp;dev.off()


#==========DENSITY LINE PLOTS================

hascounts<-output$meta_data$case[output$meta_data$WBC!=0]

ispaired <- names(which(table(hascounts)==2))

WBC = output$meta_data[output$meta_data$case %in% ispaired,]$WBC
names(WBC) = output$meta_data[output$meta_data$case %in% ispaired,]$sample_id

data.frame(
  sample_id = output$meta_data[output$meta_data$case %in% ispaired,]$sample_id,
  WBC = output$meta_data[output$meta_data$case %in% ispaired,]$WBC)

pairedprops<-props[,names(WBC)]/100

#cell # per uL
abcounts_table<-t(t(pairedprops)*WBC*1000)
abcounts <- as.data.frame.matrix(abcounts_table)




ggdf <- melt(data.frame(cluster = rownames(abcounts), abcounts, check.names = FALSE),
             id.vars = "cluster", value.name = "cells_uL", 
             variable.name = "sample_id")
ggdf$sample_id <- factor(ggdf$sample_id)
ggdf$cluster <- factor(ggdf$cluster, levels=clusterlevels)
ggdf$timepoint <- factor(output$meta_data$timepoint[match(ggdf$sample_id,output$meta_data$sample_id)], levels = timepointlevels)
ggdf$timegroup <- factor(output$meta_data$timegroup[match(ggdf$sample_id,output$meta_data$sample_id)], levels = timegrouplevels)
ggdf$response <- factor(output$meta_data$response[match(ggdf$sample_id,output$meta_data$sample_id)], levels = responselevels)
ggdf$survival <- factor(output$meta_data$survival[match(ggdf$sample_id,output$meta_data$sample_id)], levels = survivallevels)
ggdf$case <- factor(output$meta_data$case[match(ggdf$sample_id,output$meta_data$sample_id)], levels = caselevels)
ggdf$benefit <- factor(output$meta_data$benefit[match(ggdf$sample_id,output$meta_data$sample_id)], levels = benefitlevels)
ggdf$patient_id <- factor(output$meta_data$patient_id[match(ggdf$sample_id,output$meta_data$sample_id)], levels = ptlevels)
ggdf$condition <- factor(output$meta_data$condition[match(ggdf$sample_id,output$meta_data$sample_id)], levels = condlevels)
ggdf$batch <- factor(output$meta_data$batch[match(ggdf$sample_id,output$meta_data$sample_id)], levels = batchlevels)
ggdf$wbc <- factor(output$meta_data$WBC[match(ggdf$sample_id,output$meta_data$sample_id)])
ggdf$cond_tp <- paste(ggdf$condition, ggdf$timepoint, sep="_")
ggdf$cond_tp <- factor(ggdf$cond_tp, levels = c("WT_C1D1", "WT_C3D1", "Mutant_C1D1", "Mutant_C3D1"))
ggdf$resp_tp <- paste(ggdf$response, ggdf$timepoint, sep="_")
ggdf$resp_tp <- factor(ggdf$resp_tp, levels = c("PD_C1D1", "PD_C3D1", "SD_C1D1", "SD_C3D1", "PR_C1D1", "PR_C3D1"))
ggdf$surv_tp <- paste(ggdf$survival, ggdf$timepoint, sep="_")
ggdf$surv_tp <- factor(ggdf$surv_tp, levels = c("Short_C1D1", "Short_C3D1", "Long_C1D1", "Long_C3D1"))

write.csv(ggdf,"abcounts.csv")

ggp <- ggplot(ggdf, aes(x=cond_tp, y=cells_uL, group=case, col=case))+
  facet_wrap(~cluster, scales='free',ncol=6)+
  #scale_color_manual(values=c("#c93c3c","#3cc95b"))+
  scale_shape_manual(values=c(0:40))+
  geom_point()+
  geom_line(show.legend = F)+
  ggtitle("Condition (Counts)")+
  theme(axis.text.x = element_text(size=8, angle = 45, hjust=1, vjust=1, color="black"),
        axis.title.x = element_blank(),
        axis.ticks=element_line(color="black"),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=10),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(1,'lines'),
        legend.text = element_text(size=10),
        legend.key = element_rect(fill="white"))
pdf("lineplot_abcounts_all_condition+timepoint.pdf", width=11, height=8);ggp;dev.off()


+
  stat_compare_means(#comparisons = surv_tp_comparisons,
    ref.group = "Short_C1D1",
    paired=T,
    method="t.test",
    label = "p.format",
    label.x = 1.4,
    label.y.npc = "top",
    hide.ns = F,
    size=1.5))
#==========DIFFERENTIAL ABUNDANCE ANALYSIS=============


## First run a two-way ANOVA to see whether the clinical outcome interacts with the treatment effect.
## Three analyses:
## 1) proportion by treatment (before vs. after, i.e. "timepoint")
## 2) proportion by clinical outcome (stable/benefit vs. progressive/nobenefit, i.e. "condition")
## 3) whether clinical outcome has any effect on the effect of treatment on cell abundances (interaction)

aov.results<- list()
for (i in 1:length(levels(ggdf$cluster))){
  res.aov <- aov(proportion ~ timepoint*condition, data=ggdf[ggdf$cluster==levels(ggdf$cluster)[i],])
  aov.results[i] <- list(tidy(res.aov))
  names(aov.results)[i]<-as.character(ggdf$cluster[i])
}
for (i in 1:length(aov.results)){
  write.xlsx(aov.results[i],
             file=paste0(resultsDir,"/2wayANOVA.xlsx"),
             sheetName=paste(names(aov.results)[i]),
             append=T)
}

## Since clinical outcome has no bearing on the effect of treatment on the cell abundances -->
## Next, use EdgeR to compare between the two timepoints regardless of the clinical outcome

d_flowset <- output$fcs
sample_id <- factor(output$meta_data$sample_id)
#group_id <- factor(output$meta_data$condition, levels=c("NoBenefit","Benefit"))
patient_id <- factor(output$meta_data$patient_id)
timepoints <- factor(output$meta_data$timepoint)
experiment_info <- data.frame(
  patient_id,
  group_id, 
  sample_id, 
  timepoints,
  stringsAsFactors = FALSE
)
panel <- readxl::read_excel(panel_filename)
cols_markers <- c(1:length(colnames(d_flowset)))
cols_func <- match(panel$Antigen[which(panel$Functional==1)],colnames(d_flowset))
cols_lineage <- match(panel$Antigen[which(panel$Subtype==1)],colnames(d_flowset))
channel_name <- colnames(d_flowset)
marker_name <- channel_name
marker_class <- rep("none", ncol(d_flowset[[1]]))
marker_class[cols_lineage] <- "type"
marker_class <- factor(marker_class, levels = c("type", "state", "none"))
marker_info <- data.frame(
  channel_name, marker_name, marker_class, stringsAsFactors = FALSE
)
marker_info

d_se <- prepareData(d_flowset, experiment_info, marker_info)
d_se <- transformData(d_se)
d_se@elementMetadata@listData[["cluster_id"]]<-as.matrix(cell_clustering1m)
d_counts <- calcCounts(d_se)
d_counts@elementMetadata@listData[["cluster_id"]] <- factor(d_counts@NAMES)

design <- createDesignMatrix(experiment_info, cols_design = c("timepoints","patient_id"))
contrast <- createContrast(c(0, 1, rep(0,19)))
nrow(contrast) == ncol(design)
data.frame(parameters = colnames(design), contrast)#remember the first patient is omitted
res_DA <- testDA_edgeR(d_counts, design, contrast)
results <- cbind("p_val"=res_DA@elementMetadata@listData[["p_val"]],
                 "p_adj"=res_DA@elementMetadata@listData[["p_adj"]])
rownames(results)<- res_DA@elementMetadata@listData[["cluster_id"]]
write.csv(results,paste0(resultsDir,"/DAanalysis_edgeR.csv"))

#==========DIFFERENTIAL STATE PLOTS=============

exprmntbl <- data.frame(fsApply(output$fcs,exprs)[, c(output$functional_markers,"CD28")],
                        sample_id = output$sample_ids, 
                        cluster = output$cell_clustering1m) %>%
  group_by(sample_id, cluster) %>%
  summarize_all(funs(mean))
exprmntbl$cluster <- factor(exprmntbl$cluster, levels=c("Tc N", "Tc CM", "Tc EM 1", "Tc EM 2", "Tc EM 3", "Tc EM 4", "Tc EFF 1", "Tc EFF 2", "Tc 17",
                                                        "Th N 1", "Th N 2", "Th CM 1", "Th CM 2", "Th CM 3", "Th EM 1", "Th EM 2", "Th EM 3", "Th 17",
                                                        "Treg 1", "Treg 2", "DNT", "DPT", "B", "NK 1", "NK 2", "NK 3", "Non T or B"))
exprmntbl$patient_id <- factor(c(rep("001",54),rep("004",54),rep("005",54),rep("006",54),rep("007",54),rep("011",54),rep("012",54),rep("013",54),rep("014",54),rep("015",54),rep("018",54),rep("019",54),rep("022",54),rep("025",54),rep("027",54),rep("028",54),rep("030",54),rep("031",54)))
exprmntbl$timepoint <- factor(c(rep("C1D1",27),rep("C3D1",27)), levels=c("C1D1","C3D1"))
exprmntbl$response <- factor(c(rep("PD",54),rep("SD/PR",54),rep("PD",54),rep("SD/PR",54),rep("PD",162),rep("SD/PR",54),rep("PD",54),rep("SD/PR",54),rep("PD",324),rep("SD/PR",54),rep("PD",54)), levels=c("PD","SD/PR"))
exprmntbl$survival <- factor(c(rep("long",108),rep("short",54),rep("long",108),rep("short",108),rep("long",216),rep("short",378)), levels=c("short","long"))
exprmntbl$arm <- factor(c(rep("Arm B",54),rep("Arm A",54),rep("Arm B",54),rep("Arm A",216),rep("Arm B",324),rep("Arm A",54),rep("Arm B",54),rep("Arm A",108),rep("Arm B",54)), levels=c("Arm A","Arm B"))

ggdf2<-melt(exprmntbl, id.var=c("cluster","arm","sample_id","patient_id","response","survival","timepoint"))

fmlistplot <- c(output$functional_markers,"CD28");fmlistplot[10]<-"X41BB"
pdf(paste0(resultsDir,"/DSsummaryboxplotfuncmarker.pdf"),width=10,height=10)
for(i in 1:length(fmlistplot)){
  ggp <- ggplot(ggdf2[ggdf2$variable==fmlistplot[i],], aes(x=paste(survival,timepoint,sep="_"), y=value, group=patient_id, col=arm))+
    ggtitle(fmlistplot[i])+
    facet_wrap(~cluster, scales='free',ncol=6)+
    scale_color_manual(values=c("#c93c3c","#3cc95b"))+
    scale_shape_manual(values=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))+
    geom_point(aes(shape=patient_id))+
    geom_line(show.legend = F)+
    theme(axis.text.x = element_text(size=8, angle = 45, hjust=1, vjust=1, color="black"),
          axis.title.x = element_blank(),
          axis.ticks=element_line(color="black"),
          axis.line.x = element_line(size=0.5),
          axis.line.y = element_line(size=0.5),
          axis.text.y = element_text(color="black"),
          axis.title.y = element_blank(),
          strip.background = element_rect(fill=NA),
          strip.text = element_text(size=10),
          panel.background = element_rect(fill="white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.key.size = unit(1,'lines'),
          legend.text = element_text(size=10),
          legend.key = element_rect(fill="white"))
  print(ggp)
}
dev.off()

pdf(paste0(resultsDir,"/DSsummaryboxplotfuncmarker2.pdf"),width=10,height=10)
for(i in 1:length(fmlistplot)){
  ggp <- ggplot(ggdf2[ggdf2$variable==fmlistplot[i],], aes(x=paste(response,timepoint,sep="_"), y=value, group=patient_id, col=arm))+
    ggtitle(fmlistplot[i])+
    facet_wrap(~cluster, scales='free',ncol=6)+
    scale_color_manual(values=c("#c93c3c","#3cc95b"))+
    scale_shape_manual(values=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))+
    geom_point(aes(shape=patient_id))+
    geom_line(show.legend = F)+
    theme(axis.text.x = element_text(size=8, angle = 45, hjust=1, vjust=1, color="black"),
          axis.title.x = element_blank(),
          axis.ticks=element_line(color="black"),
          axis.line.x = element_line(size=0.5),
          axis.line.y = element_line(size=0.5),
          axis.text.y = element_text(color="black"),
          axis.title.y = element_blank(),
          strip.background = element_rect(fill=NA),
          strip.text = element_text(size=10),
          panel.background = element_rect(fill="white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.key.size = unit(1,'lines'),
          legend.text = element_text(size=10),
          legend.key = element_rect(fill="white"))
  print(ggp)
}
dev.off()


pdf(paste0(resultsDir,"/DSsummaryboxplotfuncmarker3.pdf"),width=10,height=12)
for(i in 1:length(fmlistplot)){
  ggp <- ggplot(ggdf2[ggdf2$variable==fmlistplot[i],], aes(x=paste(arm,timepoint,sep="_"), y=value, group=patient_id, col=response))+
    ggtitle(fmlistplot[i])+
    facet_wrap(~cluster, scales='free',ncol=6)+
    scale_color_manual(values=c("#c93c3c","#3cc95b"))+
    scale_shape_manual(values=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))+
    geom_point(aes(shape=patient_id))+
    geom_line(show.legend = F)+
    theme(axis.text.x = element_text(size=8, angle = 45, hjust=1, vjust=1, color="black"),
          axis.title.x = element_blank(),
          axis.ticks=element_line(color="black"),
          axis.line.x = element_line(size=0.5),
          axis.line.y = element_line(size=0.5),
          axis.text.y = element_text(color="black"),
          axis.title.y = element_blank(),
          strip.background = element_rect(fill=NA),
          strip.text = element_text(size=10),
          panel.background = element_rect(fill="white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.key.size = unit(1,'lines'),
          legend.text = element_text(size=10),
          legend.key = element_rect(fill="white"))
  print(ggp)
}
dev.off()



pdf(paste0(resultsDir,"/DSsummaryboxplotfuncmarker4.pdf"),width=10,height=12)
for(i in 1:length(fmlistplot)){
  ggp <- ggplot(ggdf2[ggdf2$variable==fmlistplot[i],], aes(x=paste(arm,timepoint,sep="_"), y=value, group=patient_id, col=survival))+
    ggtitle(fmlistplot[i])+
    facet_wrap(~cluster, scales='free',ncol=6)+
    scale_color_manual(values=c("#c93c3c","#3cc95b"))+
    scale_shape_manual(values=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))+
    geom_point(aes(shape=patient_id))+
    geom_line(show.legend = F)+
    theme(axis.text.x = element_text(size=8, angle = 45, hjust=1, vjust=1, color="black"),
          axis.title.x = element_blank(),
          axis.ticks=element_line(color="black"),
          axis.line.x = element_line(size=0.5),
          axis.line.y = element_line(size=0.5),
          axis.text.y = element_text(color="black"),
          axis.title.y = element_blank(),
          strip.background = element_rect(fill=NA),
          strip.text = element_text(size=10),
          panel.background = element_rect(fill="white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.key.size = unit(1,'lines'),
          legend.text = element_text(size=10),
          legend.key = element_rect(fill="white"))
  print(ggp)
}
dev.off()


#==========DIFFERENTIAL STATE ANALYSIS=============

exprmntbl_melt <- melt(exprmntbl,id.vars = c("sample_id", "cluster","condition","patient_id","timepoint"), value.name = "mean_expression", variable.name = "antigen")
exprmntbl_cluster <- dcast(exprmntbl_melt, cluster + antigen ~ sample_id, value.var = "mean_expression")
rownames(exprmntbl_cluster) <- paste0(exprmntbl_cluster$cluster, "_", exprmntbl_cluster$antigen)

## Export data table
write.csv(exprmntbl_cluster, paste0(resultsDir,'/DStable_ClinicalOutcome.csv'))

## Contrasts
k1=c(0,1,0)
k2=c(0,0,1)
contrast_names <- c("PrevsPost","StblvsProg")
K <- matrix(cbind(k1,k2), nrow = 2, byrow = T) 
rownames(K) = contrast_names

## For all paired analyses, patient-level effect will be modeled. 
## For the analysis of quantitative expression levels, batch effect will be also modeled

formula_lmer <- y ~ timepoint + condition + (1|patient_id) + (1|batch)

##Assess whether clinical outcome impacts expression
expressionda = de_wrapper(expr_mean = exprmntbl_cluster,
                          md = output[["meta_data"]][output[["meta_data"]]$sample_id,],
                          model = "lmer",
                          formula = formula_lmer,
                          K = K[2,],
                          contrast_names = contrast_names[2])
apply(expressionda$adjp < 0.1, 2, table)
write.xlsx(expressionda, paste0(resultsDir,'/DSanalysis_ClinicalOutcome.xlsx'))
write.csv(expressionda, paste0(resultsDir,'/DSanalysis_ClinicalOutcome.csv'))


##Compare Pre vs. Post Ipi/GVAX
expressionda2 = de_wrapper(expr_mean = exprmntbl_cluster,
                           md = output[["meta_data"]][output[["meta_data"]]$sample_id,],
                           model = "lmer",
                           formula = formula_lmer,
                           K = K[1,],
                           contrast_names = contrast_names[1])
apply(expressionda2$adjp < 0.1, 2, table)

## Tack on fold-changes
foldchange<-data.frame(rowMeans(exprmntbl_cluster[,seq(4,42,2)]/exprmntbl_cluster[,seq(3,42,2)]))
expressionda3<-cbind(data.frame(expressionda2),foldchange);colnames(expressionda3)[3]<-"FoldChange_PostPre"
write.xlsx(expressionda3, paste0(resultsDir,'/DSanalysis_Ipigvax.xlsx'))
write.csv(expressionda3, paste0(resultsDir,'/DSanalysis_Ipigvax.csv'))


## Plot top functional marker changes with Ipi/GVAX
filteredDS <- expressionda3[expressionda3$adjp_PrevsPost<0.1,]
filteredDS <- filteredDS[filteredDS$FoldChange_PostPre!=Inf,]#get rid of Inf (Tc17 Ki67, divide by 0 baseline)
filteredDS <- filteredDS[filteredDS$`FoldChange_PostPre`<0.8|filteredDS$`FoldChange_PostPre`>1.2,]
filteredDS[,4]<-rownames(filteredDS);colnames(filteredDS)[4]<-"celltype"
filteredDS <- filteredDS[order(filteredDS$celltype),]
filteredDS$celltype <- factor(filteredDS$celltype, levels=filteredDS$celltype)

topfuncmark <- ggplot(filteredDS, aes(x=celltype, y=FoldChange_PostPre, ))+
  geom_point(aes(size=-log2(adjp_PrevsPost)))+
  geom_segment( aes(x=celltype, xend=celltype, y=1, yend=FoldChange_PostPre))+
  scale_size(breaks=seq(1,20,4))+
  scale_y_continuous(limits = c(1,2.7),expand = c(0,0))+
  xlab("Markers within Cell Type Clusters")+
  ylab("Fold Change (Post/Pre)")+
  ggtitle("Top Functional Marker Changes with Ipi/GVAX",
          subtitle = "(FDR p value <0.05 and at Least 20% Change)")+
  theme(
    panel.background=element_blank(),
    axis.ticks.x=element_blank(),
    axis.ticks.y=element_line(size=0.25),
    axis.text=element_text(color="black"),
    axis.line=element_line(color="black", size=0.25),
    axis.text.x = element_text(size=6.8, angle = 90, hjust=1, vjust=1, color="black"),
    legend.key.size=unit(0.8,'lines'),
    legend.key = element_rect(fill="white"),
    legend.text = element_text(size=7),
    legend.title = element_text(size=7),
    legend.position = "bottom"
  )
pdf(paste0(resultsDir,"/DS_topfuncmarkers.pdf"), width=6.5, height=4);topfuncmark;dev.off()