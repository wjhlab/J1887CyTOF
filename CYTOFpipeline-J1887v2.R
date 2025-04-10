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



####CLUSTER HEATMAP FUNCTIONS ####
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

samplevels <- factor(output$meta_data$sample_id)

ptlevels <- factor(output$meta_data$patient_id)

caselevels <- factor(unique(output$meta_data$case))

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
