rm(list=ls())

devtools::load_all(".\\LouvainBayesian")
devtools::load_all(".\\monocle")
library(Seurat)
library(SPATA2)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(msigdbr)
library(cowplot)
library(clustree)
library(confuns)
library(ggsci)
library(GSVA)
library(gson)
library(GSEABase)
library(data.table)
library(spacexr)
library(tidyverse)
library(stringr)
library(patchwork)
library(STdeconvolve)
library(pals)
library(mistyR)
library(CARD)
library(distances)
library(future)
library(recipes)


###---------------------part0:load pre------------------------------------------
Load10X_Spatial_change <- function(data.dir, 
                                   filename = "filtered_feature_bc_matrix.h5", 
                                   image.name = "tissue_lowres_image.png",
                                   assay = "Spatial", slice = "slice1", 
                                   filter.matrix = TRUE, 
                                   to.upper = FALSE, image = NULL, ...) 
{
  if (length(x = data.dir) > 1) {
    warning("'Load10X_Spatial' accepts only one 'data.dir'", 
            immediate. = TRUE)
    data.dir <- data.dir[1]
  }
  if(grepl("h5", filename )){
    data <- Read10X_h5(filename = file.path(data.dir, filename), 
                       ...)
  } else {
    #data = Read10X(filtered_dir)
    data = Read10X(paste(data.dir, filename,sep = "/") )
  }
  if (to.upper) {
    rownames(x = data) <- toupper(x = rownames(x = data))
  }
  object <- CreateSeuratObject(counts = data, assay = assay)
  if (is.null(x = image)) {
    image <- Read10X_Image(image.dir = file.path(data.dir, 
                                                 "spatial"), 
                           image.name = image.name,
                           filter.matrix = filter.matrix)
  }
  else {
    if (!inherits(x = image, what = "VisiumV1")) 
      stop("Image must be an object of class 'VisiumV1'.")
  }
  image <- image[Cells(x = object)]
  DefaultAssay(object = image) <- assay
  object[[slice]] <- image
  if(image.name == "tissue_lowres_image.png") {
    object = object
  }else {
    object@images[[1]]@scale.factors$lowres = object@images[[1]]@scale.factors$hires
    }
  return(object)
}




####--------------------part1:load ST data--------------------------------------
data.dir = ".//ST_outs/"
seurat <- Load10X_Spatial_change(data.dir = data.dir,
                                 filename = "filtered_feature_bc_matrix",
                                 image.name = "tissue_hires_image.png")
seurat <- NormalizeData(seurat,normalization.method = "LogNormalize",
                        scale.factor = 10000)
seurat <- FindVariableFeatures(seurat,selection.method = "vst",nfeatures = 2000)
all.genes <- VariableFeatures(seurat)
seurat <- ScaleData(seurat,features = all.genes)
seurat <- RunPCA(seurat, features = all.genes)
seurat<-JackStraw(seurat,num.replicate=10)
seurat<-ScoreJackStraw(seurat,dims=1:20)
choose_PCA<- JackStrawPlot(seurat,dims=1:20)
choose_PCA
ndim <- 15
seurat <- FindNeighbors(seurat, dims = 1:ndim)
seq <- seq(0.1,1,by = 0.1)
for(res in seq){
  seurat <- FindClusters(seurat,resolution = res)
}
clustree(seurat, prefix = 'Spatial_snn_res.') + coord_flip()
seurat <- FindClusters(seurat, resolution = 0.6,verbose = FALSE) 
SpatialDimPlot(seurat, label.size = 3,pt.size.factor = 3,alpha = 1)
diet.seurat = Seurat::DietSeurat(seurat, graphs = "pca")
sce = as.SingleCellExperiment(diet.seurat)
colData(sce) = cbind(colData(sce), seurat@images$slice1@coordinates)
sce = spatialPreprocess(sce,platform = "ST",n.PCs = ndim,n.HVGs = 2000,
                        log.normalize = T)
q = length(unique(seurat$seurat_clusters))
sce <- spatialCluster(sce,
                      q = q, 
                      platform = "ST",
                      d = ndim, 
                      nrep = 10000, 
                      burn.in = 100,
                      init.method = "louvain",
                      model = "t",
                      gamma = 3)
seurat@meta.data = cbind(seurat@meta.data, LouvainBayesian = 
                           as.factor(sce$spatial.cluster))
Idents(seurat) <- seurat@meta.data$LouvainBayesian
color = c("#75C8AE","#DB2834","#E995C9","#FC9871","#A77852","#87CEEB","#B0DC66")
names(color) <- Idents(seurat) %>% levels()
seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:ndim)
umap_plot <- DimPlot(seurat, reduction = "umap", label = TRUE, cols = color,
                     group.by = "LouvainBayesian")
Spatial_plot<- SpatialDimPlot(seurat, "LouvainBayesian",label.size = 3, 
                              cols = color, pt.size.factor = 3,alpha = 1,
                              stroke = NA)
Spatial_plot
ggsave(
  filename = "Spatial_cluster.pdf",  
  plot = Spatial_plot,                
  width = 1120/300,                     
  height = 945/300,                      
  units = "in"                     
)
cluster_result <- data.frame(seurat$LouvainBayesian)
write.csv(cluster_result,file = "cluster_result.csv")
####--------------------part2:save color---------------------------------------
color_group <- ggplot_build(umap_plot)
colour = color_group[["data"]][[1]][["colour"]]
group =color_group[["data"]][[1]][["group"]]
info <- data.frame(cbind(colour,group))
info <- unique(arrange(info,group))
mycol <- as.character(info$colour)
mycol <- setNames(info$colour,as.character(info$group))
####--------------------part3:CNV analysis in cluster---------------------------
seurat_raw <- seurat
seurat <- subset(seurat, idents = c(1,4,5,6))
dir.create("BRCA_cnv_obj")
coords_df <- Seurat::GetTissueCoordinates(seurat) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'barcodes')
colnames(coords_df)[2:3] <- c('y','x')
spata_obj <- SPATA2::transformSeuratToSpata(seurat_object = seurat,
                                            method = 'spatial',
                                            sample_name = 'BRCA',
                                            assay_name =  'Spatial',
                                            coords_from = coords_df
)
spata_obj <- setCoordsDf(spata_obj, coords_df)
spata_obj <- runCnvAnalysis(object = spata_obj,
                            directory_cnv_folder = "BRCA_cnv_obj",
                            cnv_prefix = "Chr")
cnv_results <- getCnvResults(spata_obj)
chr_cnv <- plotCnvHeatmap(object = spata_obj, across = "LouvainBayesian", 
                          clrp_adjust = mycol)
chr_cnv
ggsave(
  filename = "chr_cnv.pdf",  
  plot = chr_cnv,                
  width = 4000/300,                     
  height = 2100/300,                      
  units = "in"                    
)
chr_cluster<- plotCnvLineplot(
  object = spata_obj,
  across = "LouvainBayesian",
  n_bins_bcsp = 1000,
  n_bins_genes = 1000,
  nrow = 3
)
chr_cluster
ggsave(
  filename = "chr_cluster.pdf",  
  plot = chr_cluster,                
  width = 4200/300,                     
  height = 2100/300,                      
  units = "in"                     
)
data_test <- as.data.frame(spata_obj@cnv[["BRCA"]][["cnv_df"]])
write.csv(data_test,file = "CNV_result.csv")
data_test2 <- data_test[,-1]
data_test2$mean <- rowMeans(data_test2)
barcodes <- as.data.frame(spata_obj@fdata[["BRCA"]][["barcodes"]])
names(barcodes) <- "barcodes"
data_test$mean <- data_test2$mean
result <- merge(barcodes,data_test,by = "barcodes",all.x = TRUE)
spata_obj@fdata[["BRCA"]][["mean"]] <- result$mean
cluster_bayes <- as.data.frame(spata_obj@fdata[["BRCA"]][["LouvainBayesian"]])
names(cluster_bayes) <- "LouvainBayesian"
data <- cbind(barcodes,cluster_bayes)
result <- merge(result,data,by = "barcodes",all.x = TRUE)
result2 <- result[,c(26,27)]
result_long <- result2 %>%
  pivot_longer(cols = c(mean), names_to = "variable", values_to = "mean")
result_long$LouvainBayesian <- factor(result_long$LouvainBayesian, 
                                      levels = c(1,4,5,6))
comparisons <- list(
  c("1", "4"),
  c("1", "5"),
  c("6", "4"),
  c("6", "5")
)
CNV_result<- ggplot(result_long, aes(x = LouvainBayesian, y = mean, 
                                     fill = LouvainBayesian)) +
  geom_violin() +
  geom_boxplot(width = 0.2 ,fill = "white",outliers.shape = NA)+
  scale_fill_manual(values=mycol)+
  scale_y_continuous(limits = c(min(result2$mean), max(result2$mean)))+
  labs(x = "", y = "CNV result") +
  theme_classic()+
  theme(legend.position = "none",
        axis.text = element_text(size=12,face="plain",color="black"))+
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif"
  )
CNV_result
ggsave(
  filename = "CNV_result.pdf",  
  plot = CNV_result,                
  width = 1010/300,                     
  height = 850/300,                      
  units = "in"                     
)

####--------------------part4:EMT score in cluster------------------------------
gene_sets <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::filter(gs_name=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
EMT <-list(EMT_signaling=unique(gene_sets$gene_symbol)) 
seurat <- AddModuleScore(seurat,features = EMT)
colnames(seurat@meta.data)[colnames(seurat@meta.data) == "Cluster1"] <- "EMT_score"
comparisons <- list(
  c("1", "4"),
  c("1", "5"),
  c("6", "4"),
  c("6", "5")
)
EMT_scores<- ggplot(seurat@meta.data, aes(x=LouvainBayesian, y=EMT_score,
                                          fill=LouvainBayesian)) +
  geom_violin() +
  geom_boxplot(width = 0.2 ,fill = "white",outliers = FALSE)+
  scale_fill_manual(values=mycol)+
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif"
  )+
  labs(x = "", y = "EMT score") +
  theme_classic()+
  theme(legend.position = "none",
        axis.text = element_text(size=12,face="plain",color="black"))
EMT_scores
ggsave(
  filename = "EMT_scores_violin.pdf",  # 输出文件名
  plot = EMT_scores,                # 要保存的图形
  width = 1010/300,                     # 图像宽度
  height = 850/300               # 分辨率，Dots Per Inch
)

####--------------------part5:metabolite pathway analysis in cluster----------------
hallmark <- read.gmt("metabolism pathway from KEGG.gmt")
hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)
matrix_all = gsva(expr = as.matrix(seurat@assays$Spatial@data), 
                  hallmark.list, 
                  kcdf="Gaussian",
                  method="ssgsea", 
                  parallel.sz=4)
matrix_all <- as.data.frame(t(matrix_all))
matrix_all <- scale(matrix_all)
LouvainBayesian <- data.frame(seurat$LouvainBayesian)
metabolite_data_result <- merge(LouvainBayesian,matrix_all,by = 0)
rownames(metabolite_data_result) <- metabolite_data_result[,1]
metabolite_data_result <- metabolite_data_result[,-1]
metabolite_data_result$seurat.LouvainBayesian <- as.factor(metabolite_data_result$seurat.LouvainBayesian)
metabolite_data_result <- metabolite_data_result %>%
  group_by(seurat.LouvainBayesian) %>%
  summarise(across(where(is.numeric), ~ mean(., na.rm = TRUE)))
metabolite_data_result <- data.frame(metabolite_data_result)
rownames(metabolite_data_result) <- paste0("cluster",
                                           metabolite_data_result$seurat.LouvainBayesian)
metabolite_data_result <- metabolite_data_result[,-1]
metabolite_data_result <- scale(metabolite_data_result)
metabolite_data_result <- data.frame(t(metabolite_data_result))
metabolite_data_result$Pathway <- rownames(metabolite_data_result)
rownames(metabolite_data_result) <- NULL
long_data <- reshape2::melt(metabolite_data_result, varnames = c("Pathway"), 
                            value.name = "Value")
p1 <- ggplot(long_data, aes(x = variable, y = Pathway, fill = Value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Heatmap of Metabolic Pathways",
       x = " ",
       y = " ",
       fill = "Value")
p1
ggsave(filename = "heatmap.pdf",
       plot = p1,
       width = 3000/300,                     
       height = 4000/300,                      
       units = "in")
data <- data.frame(seurat$LouvainBayesian)
metabolite_data_result2 <- cbind(data,matrix_all)
metabolite_data_result2 <- metabolite_data_result2[,c("seurat.LouvainBayesian",
                                                      "Sphingolipid_metabolism")]
comparisons <- list(
  c("1", "4"),
  c("1", "5"),
  c("6", "4"),
  c("6", "5")
)
Sphingolipid_metabolism <- ggplot(metabolite_data_result2, 
                                  aes(x=seurat.LouvainBayesian, 
                                      y=Sphingolipid_metabolism, 
                                      fill=seurat.LouvainBayesian)) +
  geom_boxplot(notch=TRUE)+
  scale_fill_manual(values=mycol)+
  labs(x = "", y = "Sphingolipid metabolism") +
  theme_classic()+
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif"
  )+
  theme(legend.position = "none",
        axis.text = element_text(size=12,face="plain",color="black"))
Sphingolipid_metabolism
ggsave(
  filename = "Sphingolipid metabolism.pdf",  
  plot = Sphingolipid_metabolism,                
  width = 1010/300,                     
  height = 850/300,                      
  units = "in"                     
)
write.csv(matrix_all,file = "metabolite_pathway.csv")
####--------------------part6:monocle analysis----------------------------
combin.data <- subset(seurat_raw, idents = c(1,4,5,6))
SpatialDimPlot(combin.data, label = TRUE, label.size = 3,pt.size.factor = 3)
expression_matrix = combin.data@assays$Spatial@counts
cell_metadata <- data.frame(group = combin.data[['orig.ident']],
                            clusters = Idents(combin.data))
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix), 
                              stringsAsFactors = F) 
rownames(gene_annotation) <- rownames(expression_matrix)
pd <- new("AnnotatedDataFrame", data = cell_metadata)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
HSMM <- newCellDataSet(expression_matrix,
                       phenoData = pd,
                       featureData = fd,
                       expressionFamily=negbinomial.size())
HSMM_myo <- estimateSizeFactors(HSMM)
HSMM_myo <- estimateDispersions(HSMM_myo)
diff_test_res <- differentialGeneTest(HSMM_myo,
                                      fullModelFormulaStr = '~clusters', cores = 4)
ordering_genes <- subset(diff_test_res, qval < 0.05)[,'gene_short_name']
HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)
HSMM_myo <- reduceDimension(HSMM_myo, max_components=2, method = 'DDRTree')
HSMM_myo <- orderCells(HSMM_myo)
# HSMM_myo <- orderCells(HSMM_myo, root_state = 4)
cluster_trajectory <- plot_cell_trajectory(HSMM_myo, 
                                           color_by = "clusters",cell_size =0.8) +
  scale_color_manual(values = mycol)
cluster_trajectory
ggsave(
  filename = "cluster_trajectory.pdf",  
  plot = cluster_trajectory,                
  width = 1150/300,                     
  height = 970/300,                      
  units = "in"                     
)
time_trajectory<- plot_cell_trajectory(HSMM_myo,color_by = "Pseudotime",
                                       cell_size =0.8)
time_trajectory
ggsave(
  filename = "time_trajectory.pdf",  
  plot = time_trajectory,                
  width = 1150/300,                     
  height = 970/300,                      
  units = "in"                     
)
state_trajectory<- plot_cell_trajectory(HSMM_myo,color_by = "State",cell_size =0.8)
state_trajectory
ggsave(
  filename = "state_trajectory.pdf",  
  plot = state_trajectory,                
  width = 1150/300,                     
  height = 970/300,                      
  units = "in"                     
)
plot_cell_trajectory(HSMM_myo, color_by = "clusters",cell_size =0.5) + 
  facet_wrap(~clusters, nrow = 4)
cell_Pseudotime <- data.frame(pData(HSMM_myo)$Pseudotime)
rownames(cell_Pseudotime) <- rownames(cell_metadata)
combin.data[['Pseudotime']] <- 0
combin.data[['Pseudotime']][rownames(cell_Pseudotime),] <- cell_Pseudotime
Time_spatial<- SpatialFeaturePlot(combin.data, features = c("Pseudotime"),
                                  stroke = NA,images='slice1',pt.size.factor = 3)
Time_spatial<- Time_spatial + theme(legend.position = "right")
Time_spatial
ggsave(
  filename = "Time_spatial.pdf",  
  plot = Time_spatial,                
  width = 1120/300,                     
  height = 945/300,                      
  units = "in"                    
)
branch_point_choose = 2
BEAM_res <- BEAM(HSMM_myo, branch_point = branch_point_choose,
                 progenitor_method = "duplicate",cores = 4)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
genes_branched_heatmap<- plot_genes_branched_heatmap(HSMM_myo[row.names(subset(BEAM_res,
                                                                               qval < 1e-4)),],
                                                     branch_point = branch_point_choose,
                                                     num_clusters = 4,
                                                     cores = 4,
                                                     use_gene_short_name = T,
                                                     show_rownames = T,
                                                     return_heatmap = T)
ggsave(
  filename = "genes_branched_heatmap.pdf",  
  plot = genes_branched_heatmap$ph_res,                
  width = 1150/300,                     
  height = 9700/300,                      
  units = "in"                    
)
top30_gene = BEAM_res[order(BEAM_res$qval),][1:30,'gene_short_name']
monocle_heatmap_top30 <- plot_genes_branched_heatmap(HSMM_myo[top30_gene,],
                                   branch_point = branch_point_choose,
                                   num_clusters = 4,
                                   cores = 4,
                                   use_gene_short_name = T,
                                   show_rownames = T,
                                   return_heatmap = T)
ggsave(
  filename = "monocle_heatmap_top30.pdf",  
  plot = monocle_heatmap_top30$ph_res,                
  width = 1150/300,                     
  height = 970/300,                      
  units = "in"                   
)
gene <- row.names(subset(fData(HSMM_myo),
                         gene_short_name %in% c("PSAP","CYBA")))
sample <- plot_genes_branched_pseudotime(HSMM_myo[gene,],
                               branch_point = branch_point_choose,
                               color_by = "clusters",
                               ncol = 2)+
  scale_color_manual(values = mycol)
sample
sample_data <- ggplot_build(sample)
spot_data <- sample_data$data[[1]]
names(spot_data)[1] <- "Cluster"
line_data <- sample_data$data[[2]]
unique(spot_data$Cluster)
mycol
spot_data[,1][spot_data[,1]=="#75C8AE"] <- "1"
spot_data[,1][spot_data[,1]=="#FC9871"] <- "4"
spot_data[,1][spot_data[,1]=="#A77852"] <- "5"
spot_data[,1][spot_data[,1]=="#87CEEB"] <- "6"
spot_data[,1] <- ordered(spot_data[,1], levels=c("1","4","5","6"))
spot_1 <- spot_data[spot_data$PANEL==1,]
spot_2 <- spot_data[spot_data$PANEL==2,]
line_1 <- line_data[line_data$PANEL==1,]
line_2 <- line_data[line_data$PANEL==2,]
p2 <- ggplot()+
  geom_point(data = spot_2, aes(x=x, y=y, color=Cluster), size=2, alpha=0.8)+
  theme_cowplot()+
  scale_color_manual(values = c("#75C8AE",  "#FC9871", "#A77852" ,"#87CEEB"))+
  geom_line(data = line_2[line_2$group==1,], aes(x=x, y=y,linetype = "Branch1"), lwd=1)+
  geom_line(data = line_2[line_2$group==2,], aes(x=x, y=y,linetype = "Branch2"), lwd=1)+ 
  labs(x = "Pseudotime(stretched)", y = "Expression")+
  scale_linetype_manual(name = " ",  
                        values = c(1, 3),  
                        labels = c("Cell fate 1", "Cell fate 2")) +
  ggtitle("PSAP") +  
  theme(plot.title = element_text(hjust = 0.5))  
p2
ggsave(
  filename = "monocle2_gene.pdf",  
  plot = p2,                
  width = 1150/200,                     
  height = 970/200,                      
  units = "in"                    
)
####--------------------part7:correlation PSAP with Sphingolipid metabolism--------------
PSAP_gene_expression <- log2(exprs(HSMM_myo)['PSAP',]+1)
PSAP_with_metabolite_expression <- merge(PSAP_gene_expression,matrix_all,by = 0)
df <- data.frame(PSAP = PSAP_with_metabolite_expression$x,
                 Metabolite = scale(PSAP_with_metabolite_expression$Sphingolipid_metabolism))
correlation <- ggscatter(df, x = "PSAP", 
                         y = "Metabolite",
                         add = "reg.line",
                         conf.int = TRUE,
                         cor.coef = TRUE,
                         cor.method = "pearson",
                         xlab = "PSAP expression",
                         ylab = "Sphingolipid metabolism",
                         title = "") +
  theme_classic()
correlation
ggsave(
  filename = "correlation.pdf",  
  plot = correlation,                
  width = 1150/300,                     
  height = 970/300,                      
  units = "in"                      
)
####--------------------part8:EMT score in different state----------------------
seurat <- subset(combin.data, idents = c(4,5))
seurat <- RenameIdents(seurat, `4` = "cell fate 1", `5` = "cell fate 2")
seurat$state_cluster <- Idents(seurat)
Spatial_plot<- SpatialDimPlot(seurat, label.size = 3, 
                              cols = color, pt.size.factor = 3,alpha = 1,
                              stroke = NA)
Spatial_plot
gene_sets <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::filter(gs_name=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
EMT <-list(EMT_signaling=unique(gene_sets$gene_symbol)) 
seurat <- AddModuleScore(seurat,features = EMT)
colnames(seurat@meta.data)[colnames(seurat@meta.data) == "Cluster1"] <- "EMT_score"
comparisons <- list(
  c("cell fate 1", "cell fate 2")
)

EMT_scores<- ggplot(seurat@meta.data,
                    aes(x=state_cluster, y=EMT_score,fill=state_cluster)) +
  geom_violin() +
  geom_boxplot(width = 0.2 ,fill = "white",outliers = FALSE)+
  scale_fill_manual(values=c("#D9A896","#A77852"))+
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif"
  )+
  labs(x = "", y = "EMT score") +
  theme_classic()+
  theme(legend.position = "none",
        axis.text = element_text(size=12,face="plain",color="black"))
EMT_scores
ggsave(
  filename = "EMT_scores_violin.pdf",  # 输出文件名
  plot = EMT_scores,                # 要保存的图形
  width = 800/300,                     # 图像宽度
  height = 850/300               # 分辨率，Dots Per Inch
)






####--------------------part9:different state colour-----------------------
state_cluster <- data.frame(HSMM_myo@phenoData@data)
cluster_LouvainBayesian <- data.frame(combin.data$LouvainBayesian)
combin.data@meta.data = cbind(combin.data@meta.data, 
                              state_cluster = as.factor(state_cluster$State))
Idents(combin.data) <- combin.data@meta.data$state_cluster
default_colors <- scales::hue_pal()(length(unique(combin.data@meta.data$state_cluster)))
print(default_colors)
####--------------------part10:metabolite pathways in different state------------------
metabolite_data_result_state <- merge(state_cluster,matrix_all,by = 0)
rownames(metabolite_data_result_state) <- metabolite_data_result_state[,1]
metabolite_data_result_state <- metabolite_data_result_state[,-c(1:5)]
metabolite_data_result_state$State <- as.factor(metabolite_data_result_state$State)
metabolite_data_result_state <- metabolite_data_result_state %>%
  group_by(State) %>%
  summarise(across(where(is.numeric), ~ mean(., na.rm = TRUE)))
metabolite_data_result_state <- data.frame(metabolite_data_result_state)
rownames(metabolite_data_result_state) <- paste0("state",
                                                 metabolite_data_result_state$State)
metabolite_data_result_state <- metabolite_data_result_state[,-1]
metabolite_data_result_state <- scale(metabolite_data_result_state)
metabolite_data_result_state <- data.frame(t(metabolite_data_result_state))
metabolite_data_result_state$Pathway <- rownames(metabolite_data_result_state)
rownames(metabolite_data_result_state) <- NULL
long_data_state <- melt(metabolite_data_result_state, varnames = c("Pathway"), 
                        value.name = "Value")
p1 <- ggplot(long_data_state, aes(x = variable, y = Pathway, fill = Value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Heatmap of Metabolic Pathways",
       x = " ",
       y = " ",
       fill = "Value")
p1
ggsave(filename = "heatmap2.pdf",
       plot = p1,
       width = 2500/300,                     
       height = 3000/300,                      
       units = "in")
state_cluster <- data.frame(combin.data$state_cluster)
metabolite_data_result_state <- merge(state_cluster,matrix_all,by = 0)
rownames(metabolite_data_result_state) <- metabolite_data_result_state[,1]
metabolite_data_result_state <- metabolite_data_result_state[,-1]
metabolite_data_result_state$combin.data.state_cluster <- as.factor(metabolite_data_result_state$combin.data.state_cluster)
metabolite_data_result_state2 <- metabolite_data_result_state[,c("combin.data.state_cluster","Sphingolipid_metabolism")]
Sphingolipid_metabolism_state <- ggplot(metabolite_data_result_state2, 
                                        aes(x=combin.data.state_cluster, 
                                            y=Sphingolipid_metabolism,
                                            fill=combin.data.state_cluster)) +
  geom_boxplot(notch=TRUE)+
  labs(x = "", y = "Sphingolipid metabolism") +
  theme_classic()+
  theme(legend.position = "none",
        axis.text = element_text(size=12,face="plain",color="black"))
Sphingolipid_metabolism_state
ggsave(
  filename = "Sphingolipid metabolism2.pdf",  
  plot = Sphingolipid_metabolism_state,                
  width = 1010/300,                     
  height = 850/300,                      
  units = "in"                    
)
####--------------------part11:Sphingolipid metabolism in state4 and state5---------------
state_cluster2 <- data.frame(combin.data$state_cluster)
state_cluster2 <- state_cluster2 %>%
  filter(combin.data$state_cluster %in% c(4, 5))
metabolite_data_result_state_choose <- merge(state_cluster2,matrix_all,by = 0)
rownames(metabolite_data_result_state_choose) <- metabolite_data_result_state_choose[,1]
metabolite_data_result_state_choose <- metabolite_data_result_state_choose[,-1]
metabolite_data_result_state_choose$combin.data.state_cluster <- as.factor(metabolite_data_result_state_choose$combin.data.state_cluster)
metabolite_data_result_state_choose <- metabolite_data_result_state_choose[,c("combin.data.state_cluster","Sphingolipid_metabolism")]
metabolite_data_result_state_choose$combin.data.state_cluster <- case_when(
  metabolite_data_result_state_choose$combin.data.state_cluster == 4 ~ "cell fate 1",
  metabolite_data_result_state_choose$combin.data.state_cluster == 5 ~ "cell fate 2"
)
comparisons <- list(
  c("cell fate 1", "cell fate 2")
)
Sphingolipid_metabolism_choose_state <- ggplot(metabolite_data_result_state_choose, 
                                               aes(x=combin.data.state_cluster, 
                                                   y=Sphingolipid_metabolism,
                                                   fill=combin.data.state_cluster)) +
  geom_boxplot(notch=TRUE)+
  scale_fill_manual(values=c("#D9A896","#A77852"))+
  labs(x = "", y = "Sphingolipid metabolism") +
  theme_classic()+
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif"
  )+
  theme(legend.position = "none",
        axis.text = element_text(size=12,face="plain",color="black"))
Sphingolipid_metabolism_choose_state
ggsave(
  filename = "Sphingolipid metabolism3.pdf",  
  plot = Sphingolipid_metabolism_choose_state,                
  width = 800/300,                     
  height = 850/300,                      
  units = "in"                    
)
####--------------------part12:RCTD Analysis ------------------------------
seurat_data <- Read10X(data.dir = "./GSE167036")
data_seurat <- CreateSeuratObject(counts = seurat_data,
                                  min.features = 200,
                                  min.cells = 3)
scRNA <- subset(data_seurat,idents = c("CA1","CA2","CA3","CA4","CA5","CA6",
                                       "CA7","CA8"))
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 3000)
all.genes <- VariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, features = all.genes)
scRNA <- RunPCA(scRNA, features = all.genes)
cell_type <- read.csv(file = "GSE167036_meta.csv",row.names = 1)
cell_type_select <- cell_type[cell_type$orig.ident %in% c("CA1", "CA2","CA3",
                                                          "CA4","CA5","CA6",
                                                          "CA7","CA8"), ]
scRNA@meta.data = cbind(scRNA@meta.data, cell_type = as.factor(cell_type_select$main_label))
Idents(scRNA) <- scRNA@meta.data$cell_type
Idents(scRNA) <- "celltype"
counts <- scRNA[["RNA"]]$counts
cluster <- as.factor(scRNA$cell_type)
names(cluster) <- colnames(scRNA)
nUMI <- scRNA$nCount_RNA
names(nUMI) <- colnames(scRNA)
reference <- Reference(counts, cluster, nUMI)
counts <- seurat_raw[["Spatial"]]$counts
coords <- GetTissueCoordinates(seurat_raw)
colnames(coords) <- c("x", "y")
coords[is.na(colnames(coords))] <- NULL
query <- SpatialRNA(coords, counts, colSums(counts))
RCTD <- create.RCTD(query, reference, max_cores = 4)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")
barcodes <- colnames(RCTD@spatialRNA@counts)
weights <- RCTD@results$weights
norm_weights <- normalize_weights(weights)
RCTD_color <- c("#FF79DE","#FE634F","#00A3FD","#F89C2B","#02CF97","#CA2979",
                "#76B1FD","#DBAD27","#FE8984")
plt <- vizAllTopics(theta = as.matrix(norm_weights),
                    pos = coords,
                    topicOrder=seq(ncol(norm_weights)),
                    topicCols=RCTD_color,
                    groups = NA,
                    group_cols = NA,
                    r = 17, 
                    lwd = 0.01,
                    showLegend = TRUE,
                    plotTitle = "")
plt <- plt + coord_flip() + scale_x_reverse()
plt
ggsave("RCTD_result.pdf",
       plot = plt,
       width = 8,                     
       height = 6*1.07,
       bg = "white" )
expression_matrix <- data.frame(as.matrix(norm_weights))
colnames(expression_matrix) <- c("B", "CAFs","CD4 T", "CD8 T", "Epithelial Cells",
                                 "KI67 Cells", "Myeloid Cells","NK","TECs" )
cluster_result <- data.frame(seurat$LouvainBayesian)
result_cluster_percent <- merge(cluster_result,expression_matrix,by = 0)
rownames(result_cluster_percent) <- result_cluster_percent$Row.names
result_cluster_percent <- result_cluster_percent[,-1]
result_cluster_percent$seurat.LouvainBayesian <- as.factor(result_cluster_percent$seurat.LouvainBayesian)
result_cluster_percent_mean_values <- result_cluster_percent %>%
  group_by(seurat.LouvainBayesian) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))
result_cluster_percent_mean_values <- result_cluster_percent_mean_values %>%
  pivot_longer(
    cols = -seurat.LouvainBayesian,  
    names_to = "cell_type",  
    values_to = "percent"  
  )
result_cluster_percent_mean_values <- result_cluster_percent_mean_values %>%  
  group_by(seurat.LouvainBayesian) %>%  
  mutate(Percentage = percent / sum(percent) * 100)
write.csv(expression_matrix,file = "RCTD_result.csv")
p1_cluster_percent<- ggplot(data=result_cluster_percent_mean_values,
                            aes(x=seurat.LouvainBayesian,
                                y=Percentage,
                                fill=cell_type))+
  geom_bar(stat = "identity",
           width = 0.8,   
           color="black", 
           linewidth=0.3)+ 
  theme_bw()+
  scale_fill_manual(values=RCTD_color)+
  scale_y_continuous(expand = expansion(mult=c(0,0.05),add=c(0,0)))+   
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5))+
  labs(y = "Percentage(%)")
p1_cluster_percent
ggsave(
  filename = "cell_type_percent_in_cluster.pdf",  
  plot = p1_cluster_percent,                
  width = 1010/200,                     
  height = 850/200,                      
  units = "in"                   
)
state_result <- data.frame(combin.data$state_cluster)
result_state_percent <- merge(state_result,expression_matrix,by = 0)
rownames(result_state_percent) <- result_state_percent$Row.names
result_state_percent <- result_state_percent[,-1]
result_state_percent$combin.data.state_cluster <- as.factor(result_state_percent$combin.data.state_cluster)
result_state_percent_mean_values <- result_state_percent %>%
  group_by(combin.data.state_cluster) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))
result_state_percent_mean_values <- result_state_percent_mean_values %>%
  pivot_longer(
    cols = -combin.data.state_cluster,  
    names_to = "cell_type",  
    values_to = "percent"  
  )
result_state_percent_mean_values <- result_state_percent_mean_values %>%  
  group_by(combin.data.state_cluster) %>%  
  mutate(percent = percent / sum(percent) * 100)
p1_state_percent<- ggplot(data=result_state_percent_mean_values,
                          aes(x=combin.data.state_cluster,
                              y=percent,
                              fill=cell_type))+
  geom_bar(stat = "identity",
           width = 0.8,   
           color="black", 
           linewidth=0.3)+ 
  theme_bw()+
  scale_fill_manual(values=RCTD_color)+
  scale_y_continuous(expand = expansion(mult=c(0,0.05),add=c(0,0)))+   
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5))+
  labs(y = "Percentage(%)")
p1_state_percent
ggsave(
  filename = "cell_type_percent_in_state.pdf",  
  plot = p1_state_percent,                
  width = 1010/200,                     
  height = 850/200,                      
  units = "in"                    
)

##########The proportion of cell types in different state clusters
data_percent_cluster <- result_state_percent

####B cells
data_percent_cluster_B <- data_percent_cluster[,c("combin.data.state_cluster","B")]
data_percent_cluster_B <- data_percent_cluster_B %>%
  pivot_longer(cols = c(B), names_to = "variable",values_to = "B")
colnames(data_percent_cluster_B) <- c("state_cluster","variable","mean")
data_percent_cluster_B <- data_percent_cluster_B %>%
  filter(state_cluster %in% c("4", "5"))
data_percent_cluster_B$state_cluster <- factor(data_percent_cluster_B$state_cluster, 
                                               levels = c("4","5"))
data_percent_cluster_B$state_cluster <- case_when(
  data_percent_cluster_B$state_cluster == 4 ~ "cell fate 1",
  data_percent_cluster_B$state_cluster == 5 ~ "cell fate 2"
)
comparisons <- list(
  c("cell fate 1", "cell fate 2")
)
plt_data_percent_cluster_B <- ggplot(data_percent_cluster_B, 
                                     aes(x = state_cluster, 
                                         y = mean, 
                                         fill = state_cluster)) +
  geom_boxplot(notch=TRUE)+
  scale_fill_manual(values=c("#D9A896","#A77852"))+
  labs(x = "", y = "Percentage(%)") +
  theme_classic()+
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif"
  )+
  theme(legend.position = "none",
        axis.text = element_text(size=12,face="plain",color="black"))
plt_data_percent_cluster_B
ggsave(
  filename = "data_percent_cluster_B.pdf",  
  plot = plt_data_percent_cluster_B,                
  width = 800/300,                     
  height = 850/300,                     
  units = "in"                     
)

####CAFs 
data_percent_cluster_CAFs <- data_percent_cluster[,c("combin.data.state_cluster","CAFs")]
data_percent_cluster_CAFs <- data_percent_cluster_CAFs %>%
  pivot_longer(cols = c(CAFs), names_to = "variable",values_to = "CAFs")
colnames(data_percent_cluster_CAFs) <- c("state_cluster","variable","mean")
data_percent_cluster_CAFs <- data_percent_cluster_CAFs %>%
  filter(state_cluster %in% c("4", "5"))
data_percent_cluster_CAFs$state_cluster <- factor(data_percent_cluster_CAFs$state_cluster, 
                                                  levels = c("4","5"))
data_percent_cluster_CAFs$state_cluster <- case_when(
  data_percent_cluster_CAFs$state_cluster == 4 ~ "cell fate 1",
  data_percent_cluster_CAFs$state_cluster == 5 ~ "cell fate 2"
)
comparisons <- list(
  c("cell fate 1", "cell fate 2")
)
plt_data_percent_cluster_CAFs <- ggplot(data_percent_cluster_CAFs, 
                                        aes(x = state_cluster, 
                                            y = mean, 
                                            fill = state_cluster)) +
  geom_boxplot(notch=TRUE)+
  scale_fill_manual(values=c("#D9A896","#A77852"))+
  labs(x = "", y = "Percentage(%)") +
  theme_classic()+
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif"
  )+
  theme(legend.position = "none",
        axis.text = element_text(size=12,face="plain",color="black"))
plt_data_percent_cluster_CAFs
ggsave(
  filename = "data_percent_cluster_CAFs.pdf",  
  plot = plt_data_percent_cluster_CAFs,                
  width = 800/300,                     
  height = 850/300,                
  units = "in"                     
)

####CD4 T 
data_percent_cluster_CD4T <- data_percent_cluster[,c("combin.data.state_cluster","CD4 T")]
data_percent_cluster_CD4T <- data_percent_cluster_CD4T %>%
  pivot_longer(cols = c('CD4 T'), names_to = "variable",values_to = "CD4 T")
colnames(data_percent_cluster_CD4T) <- c("state_cluster","variable","mean")
data_percent_cluster_CD4T <- data_percent_cluster_CD4T %>%
  filter(state_cluster %in% c("4", "5"))
data_percent_cluster_CD4T$state_cluster <- factor(data_percent_cluster_CD4T$state_cluster,
                                                  levels = c("4","5"))
data_percent_cluster_CD4T$state_cluster <- case_when(
  data_percent_cluster_CD4T$state_cluster == 4 ~ "cell fate 1",
  data_percent_cluster_CD4T$state_cluster == 5 ~ "cell fate 2"
)
comparisons <- list(
  c("cell fate 1", "cell fate 2")
)
plt_data_percent_cluster_CD4T <- ggplot(data_percent_cluster_CD4T, 
                                        aes(x = state_cluster, 
                                            y = mean, 
                                            fill = state_cluster)) +
  geom_boxplot(notch=TRUE)+
  scale_fill_manual(values=c("#D9A896","#A77852"))+
  labs(x = "", y = "Percentage(%)") +
  theme_classic()+
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif"
  )+
  theme(legend.position = "none",
        axis.text = element_text(size=12,face="plain",color="black"))
plt_data_percent_cluster_CD4T
ggsave(
  filename = "data_percent_cluster_CD4T.pdf",  
  plot = plt_data_percent_cluster_CD4T,                
  width = 800/300,                     
  height = 850/300,            
  units = "in"                     
)
####CD8 T 
data_percent_cluster_CD8T <- data_percent_cluster[,c("combin.data.state_cluster","CD8 T")]
data_percent_cluster_CD8T <- data_percent_cluster_CD8T %>%
  pivot_longer(cols = c('CD8 T'), names_to = "variable",values_to = "CD8 T")
colnames(data_percent_cluster_CD8T) <- c("state_cluster","variable","mean")
data_percent_cluster_CD8T <- data_percent_cluster_CD8T %>%
  filter(state_cluster %in% c("4", "5"))
data_percent_cluster_CD8T$state_cluster <- factor(data_percent_cluster_CD8T$state_cluster,
                                                  levels = c("4","5"))
data_percent_cluster_CD8T$state_cluster <- case_when(
  data_percent_cluster_CD8T$state_cluster == 4 ~ "cell fate 1",
  data_percent_cluster_CD8T$state_cluster == 5 ~ "cell fate 2"
)
comparisons <- list(
  c("cell fate 1", "cell fate 2")
)
plt_data_percent_cluster_CD8T <- ggplot(data_percent_cluster_CD8T, 
                                        aes(x = state_cluster, 
                                            y = mean, 
                                            fill = state_cluster)) +
  geom_boxplot(notch=TRUE)+
  scale_fill_manual(values=c("#D9A896","#A77852"))+
  labs(x = "", y = "Percentage(%)") +
  theme_classic()+
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif"
  )+
  theme(legend.position = "none",
        axis.text = element_text(size=12,face="plain",color="black"))
plt_data_percent_cluster_CD8T
ggsave(
  filename = "data_percent_cluster_CD8T.pdf",  
  plot = plt_data_percent_cluster_CD8T,                
  width = 800/300,                     
  height = 850/300,
  units = "in"                     
)
####Epithelial Cells 
data_percent_cluster_Epithelial_Cells <- data_percent_cluster[,c("combin.data.state_cluster","Epithelial Cells")]
data_percent_cluster_Epithelial_Cells <- data_percent_cluster_Epithelial_Cells %>%
  pivot_longer(cols = c('Epithelial Cells'), names_to = "variable",
               values_to = "Epithelial Cells")
colnames(data_percent_cluster_Epithelial_Cells) <- c("state_cluster","variable","mean")
data_percent_cluster_Epithelial_Cells <- data_percent_cluster_Epithelial_Cells %>%
  filter(state_cluster %in% c("4", "5"))
data_percent_cluster_Epithelial_Cells$state_cluster <- factor(data_percent_cluster_Epithelial_Cells$state_cluster, levels = c("4","5"))
data_percent_cluster_Epithelial_Cells$state_cluster <- case_when(
  data_percent_cluster_Epithelial_Cells$state_cluster == 4 ~ "cell fate 1",
  data_percent_cluster_Epithelial_Cells$state_cluster == 5 ~ "cell fate 2"
)
comparisons <- list(
  c("cell fate 1", "cell fate 2")
)
plt_data_percent_cluster_Epithelial_Cells <- ggplot(data_percent_cluster_Epithelial_Cells,
                                                    aes(x = state_cluster, 
                                                        y = mean, 
                                                        fill = state_cluster)) +
  geom_boxplot(notch=TRUE)+
  scale_fill_manual(values=c("#D9A896","#A77852"))+
  labs(x = "", y = "Percentage(%)") +
  theme_classic()+
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif"
  )+
  theme(legend.position = "none",
        axis.text = element_text(size=12,face="plain",color="black"))
plt_data_percent_cluster_Epithelial_Cells
ggsave(
  filename = "data_percent_cluster_Epithelial_Cells.pdf",  
  plot = plt_data_percent_cluster_Epithelial_Cells,                
  width = 800/300,                     
  height = 850/300,            
  units = "in"                     
)
####KI67 Cells
data_percent_cluster_KI67_Cells <- data_percent_cluster[,c("combin.data.state_cluster",
                                                           "KI67 Cells")]
data_percent_cluster_KI67_Cells <- data_percent_cluster_KI67_Cells %>%
  pivot_longer(cols = c('KI67 Cells'), names_to = "variable",values_to = "KI67 Cells")
colnames(data_percent_cluster_KI67_Cells) <- c("state_cluster","variable","mean")
data_percent_cluster_KI67_Cells <- data_percent_cluster_KI67_Cells %>%
  filter(state_cluster %in% c("4", "5"))
data_percent_cluster_KI67_Cells$state_cluster <- factor(data_percent_cluster_KI67_Cells$state_cluster, levels = c("4","5"))
data_percent_cluster_KI67_Cells$state_cluster <- case_when(
  data_percent_cluster_KI67_Cells$state_cluster == 4 ~ "cell fate 1",
  data_percent_cluster_KI67_Cells$state_cluster == 5 ~ "cell fate 2"
)
comparisons <- list(
  c("cell fate 1", "cell fate 2")
)
plt_data_percent_cluster_KI67_Cells <- ggplot(data_percent_cluster_KI67_Cells, 
                                              aes(x = state_cluster, 
                                                  y = mean, 
                                                  fill = state_cluster)) +
  geom_boxplot(notch=TRUE)+
  scale_fill_manual(values=c("#D9A896","#A77852"))+
  labs(x = "", y = "Percentage(%)") +
  theme_classic()+
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif"
  )+
  theme(legend.position = "none",
        axis.text = element_text(size=12,face="plain",color="black"))
plt_data_percent_cluster_KI67_Cells
ggsave(
  filename = "data_percent_cluster_KI67_Cells.pdf",  
  plot = plt_data_percent_cluster_KI67_Cells,                
  width = 800/300,                     
  height = 850/300,          
  units = "in"                     
)
####Myeloid Cells
data_percent_cluster_Myeloid_Cells <- data_percent_cluster[,c("combin.data.state_cluster",
                                                              "Myeloid Cells")]
data_percent_cluster_Myeloid_Cells <- data_percent_cluster_Myeloid_Cells %>%
  pivot_longer(cols = c('Myeloid Cells'), names_to = "variable",values_to = "Myeloid Cells")
colnames(data_percent_cluster_Myeloid_Cells) <- c("state_cluster","variable","mean")
data_percent_cluster_Myeloid_Cells <- data_percent_cluster_Myeloid_Cells %>%
  filter(state_cluster %in% c("4", "5"))
data_percent_cluster_Myeloid_Cells$state_cluster <- factor(data_percent_cluster_Myeloid_Cells$state_cluster, levels = c("4","5"))
data_percent_cluster_Myeloid_Cells$state_cluster <- case_when(
  data_percent_cluster_Myeloid_Cells$state_cluster == 4 ~ "cell fate 1",
  data_percent_cluster_Myeloid_Cells$state_cluster == 5 ~ "cell fate 2"
)
comparisons <- list(
  c("cell fate 1", "cell fate 2")
)
plt_data_percent_cluster_Myeloid_Cells <- ggplot(data_percent_cluster_Myeloid_Cells, 
                                                 aes(x = state_cluster, 
                                                     y = mean, 
                                                     fill = state_cluster)) +
  geom_boxplot(notch=TRUE)+
  scale_fill_manual(values=c("#D9A896","#A77852"))+
  labs(x = "", y = "Percentage(%)") +
  theme_classic()+
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif"
  )+
  theme(legend.position = "none",
        axis.text = element_text(size=12,face="plain",color="black"))
plt_data_percent_cluster_Myeloid_Cells
ggsave(
  filename = "data_percent_cluster_Myeloid_Cells.pdf",  
  plot = plt_data_percent_cluster_Myeloid_Cells,                
  width = 800/300,                     
  height = 850/300,             
  units = "in"                     
)
####NK
data_percent_cluster_NK <- data_percent_cluster[,c("combin.data.state_cluster","NK")]
data_percent_cluster_NK <- data_percent_cluster_NK %>%
  pivot_longer(cols = c('NK'), names_to = "variable",values_to = "NK")
colnames(data_percent_cluster_NK) <- c("state_cluster","variable","mean")
data_percent_cluster_NK <- data_percent_cluster_NK %>%
  filter(state_cluster %in% c("4", "5"))
data_percent_cluster_NK$state_cluster <- factor(data_percent_cluster_NK$state_cluster, 
                                                levels = c("4","5"))
data_percent_cluster_NK$state_cluster <- case_when(
  data_percent_cluster_NK$state_cluster == 4 ~ "cell fate 1",
  data_percent_cluster_NK$state_cluster == 5 ~ "cell fate 2"
)
comparisons <- list(
  c("cell fate 1", "cell fate 2")
)
plt_data_percent_cluster_NK <- ggplot(data_percent_cluster_NK, 
                                      aes(x = state_cluster,
                                          y = mean,
                                          fill = state_cluster)) +
  geom_boxplot(notch=TRUE)+
  scale_fill_manual(values=c("#D9A896","#A77852"))+
  labs(x = "", y = "Percentage(%)") +
  theme_classic()+
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif"
  )+
  theme(legend.position = "none",
        axis.text = element_text(size=12,face="plain",color="black"))
plt_data_percent_cluster_NK
ggsave(
  filename = "data_percent_cluster_NK.pdf",  
  plot = plt_data_percent_cluster_NK,                
  width = 800/300,                     
  height = 850/300,              
  units = "in"                     
)
####TECs
data_percent_cluster_TECs <- data_percent_cluster[,c("combin.data.state_cluster","TECs")]
data_percent_cluster_TECs <- data_percent_cluster_TECs %>%
  pivot_longer(cols = c('TECs'), names_to = "variable",values_to = "TECs")
colnames(data_percent_cluster_TECs) <- c("state_cluster","variable","mean")
data_percent_cluster_TECs <- data_percent_cluster_TECs %>%
  filter(state_cluster %in% c("4", "5"))
data_percent_cluster_TECs$state_cluster <- factor(data_percent_cluster_TECs$state_cluster, 
                                                  levels = c("4","5"))
data_percent_cluster_TECs$state_cluster <- case_when(
  data_percent_cluster_TECs$state_cluster == 4 ~ "cell fate 1",
  data_percent_cluster_TECs$state_cluster == 5 ~ "cell fate 2"
)
comparisons <- list(
  c("cell fate 1", "cell fate 2")
)
plt_data_percent_cluster_TECs <- ggplot(data_percent_cluster_TECs, 
                                        aes(x = state_cluster, 
                                            y = mean, 
                                            fill = state_cluster)) +
  geom_boxplot(notch=TRUE)+
  scale_fill_manual(values=c("#D9A896","#A77852"))+
  labs(x = "", y = "Percentage(%)") +
  theme_classic()+
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif"
  )+
  theme(legend.position = "none",
        axis.text = element_text(size=12,face="plain",color="black"))
plt_data_percent_cluster_TECs
ggsave(
  filename = "data_percent_cluster_TECs.pdf",  
  plot = plt_data_percent_cluster_TECs,                
  width = 800/300,                     
  height = 850/300,              
  units = "in"                     
)

####--------------------part13::MISTy-------------------------------------------
composition = as.matrix(norm_weights)
colnames(composition) <- c("Bcell","CAFs","CD4T","CD8T","EpithelialCells","KI67cells",
                           "Myeloidcells","NK","TECs")
###MISTy in cell fate 1
seurat <- combin.data
geometry <- GetTissueCoordinates(seurat, cols = c("imagerow", "imagecol"), scale = NULL)
data_cluster <- data.frame(seurat$state_cluster)
composition2 <- merge(composition,data_cluster,by = 0)
geometry2 <- merge(geometry,data_cluster,by = 0)
composition2 <- composition2[composition2$seurat.state_cluster %in% c("4"),]
geometry2 <- geometry2[geometry2$seurat.state_cluster %in% c("4"),]
rownames(composition2) <- composition2$Row.names
rownames(geometry2) <- geometry2$Row.names
composition <- composition2[,-c(1,11)]
geometry <- geometry2[,-c(1,4)]
geom_dist <- as.matrix(distances(geometry))
dist_nn <- apply(geom_dist, 1, function(x) (sort(x)[2]))
paraview_radius <- ceiling(mean(dist_nn+ sd(dist_nn)))
GBM_views <- create_initial_view(composition) %>%  
  add_paraview(geometry, l= paraview_radius, family = "gaussian") 
dir.create("./cellfate1", recursive = TRUE, showWarnings = FALSE)
run_misty(GBM_views, "./cellfate1")
misty_results <- collect_results("./cellfate1")
misty_results %>%  
  plot_improvement_stats("multi.R2") %>%
  plot_improvement_stats("gain.R2")
misty_results %>% 
  plot_interaction_heatmap(view = "intra", clean = F)
misty_results$importances.aggregated %>%  
  filter(view == "intra", Predictor == "CAFs") %>%
  arrange(-Importance)
paraview_radius
misty_results %>% plot_interaction_heatmap(view = "para.194", clean = F,
                                           trim = 0.05, trim.measure = "gain.R2",
                                           cutoff = 0.5)

pdf("cell_fate_1.pdf", width = 7, height = 7) 
par(mar = c(5, 5, 4, 2)) 
misty_results %>% plot_interaction_communities("para.194")
box(col = "black", lwd = 2)  
dev.off()

###MISTy in cell fate 2
seurat <- combin.data
geometry <- GetTissueCoordinates(seurat, cols = c("imagerow", "imagecol"), scale = NULL)
data_cluster <- data.frame(seurat$state_cluster)
composition2 <- merge(composition,data_cluster,by = 0)
geometry2 <- merge(geometry,data_cluster,by = 0)
composition2 <- composition2[composition2$seurat.state_cluster %in% c("5"),]
geometry2 <- geometry2[geometry2$seurat.state_cluster %in% c("5"),]
rownames(composition2) <- composition2$Row.names
rownames(geometry2) <- geometry2$Row.names
composition <- composition2[,-c(1,11)]
geometry <- geometry2[,-c(1,4)]
geom_dist <- as.matrix(distances(geometry))
dist_nn <- apply(geom_dist, 1, function(x) (sort(x)[2]))
paraview_radius <- ceiling(mean(dist_nn+ sd(dist_nn)))
GBM_views <- create_initial_view(composition) %>%  
  add_paraview(geometry, l= paraview_radius, family = "gaussian") 
dir.create("./cellfate2", recursive = TRUE, showWarnings = FALSE)
run_misty(GBM_views, "./cellfate2")
misty_results <- collect_results("./cellfate2")
misty_results %>%  
  plot_improvement_stats("multi.R2") %>%
  plot_improvement_stats("gain.R2")
misty_results %>% 
  plot_interaction_heatmap(view = "intra", clean = F)
misty_results$importances.aggregated %>%  
  filter(view == "intra", Predictor == "CAFs") %>%
  arrange(-Importance)
paraview_radius
misty_results %>% plot_interaction_heatmap(view = "para.194", clean = F,
                                           trim = 0.05, trim.measure = "gain.R2",
                                           cutoff = 0.5)

pdf("cell_fate_2.pdf", width = 7, height = 7) 
par(mar = c(5, 5, 4, 2)) 
misty_results %>% plot_interaction_communities("para.194")
box(col = "black", lwd = 2)  
dev.off()


####--------------------part13:save result -------------------------------------
save.image(file = "ST.RData")
