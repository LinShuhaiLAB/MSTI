rm(list=ls())

devtools::load_all(".\\LouvainBayesian")
library(Seurat)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(msigdbr)
library(cowplot)
library(clustree)
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


####--------------------part1:load SM data--------------------------------------
data.dir = ".//SM_outs/"
seurat <- Load10X_Spatial_change(data.dir = data.dir,
                                 filename = "filtered_feature_bc_matrix",
                                 image.name = "tissue_lowres_image.png")
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize",
                        scale.factor = 10000)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 350)
all.genes <- VariableFeatures(seurat)
seurat <- ScaleData(seurat, features = all.genes)
seurat <- RunPCA(seurat, features = all.genes)
seurat<-JackStraw(seurat, num.replicate=10)
seurat<-ScoreJackStraw(seurat, dims=1:20)
JackStrawPlot(seurat, dims=1:20)
ndim <- 9
seurat <- FindNeighbors(seurat, dims = 1:ndim)
seq <- seq(0.1,1, by = 0.1)
for(res in seq){
  seurat <- FindClusters(seurat, resolution = res)
}
clustree(seurat, prefix = 'Spatial_snn_res.') + coord_flip()
seurat <- FindClusters(seurat, resolution = 0.4, verbose = FALSE)
Spatial_plot <- SpatialDimPlot(seurat, label.size = 3, pt.size.factor = 1,
                               alpha = 1, stroke = NA)
diet.seurat = Seurat::DietSeurat(seurat, graphs = "pca")
sce = as.SingleCellExperiment(diet.seurat)
colData(sce) = cbind(colData(sce), seurat@images$slice1@coordinates)
sce = spatialPreprocess(sce,platform = "ST", n.PCs = ndim, n.HVGs = 2000, log.normalize = T)
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
seurat@meta.data = cbind(seurat@meta.data, LouvainBayesian = as.factor(sce$spatial.cluster))
Idents(seurat) <- seurat@meta.data$LouvainBayesian
seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:ndim)
umap_plot <- DimPlot(seurat, reduction = "umap", label = TRUE,
                     group.by = "LouvainBayesian")
umap_plot
ggsave(
  filename = "UMAP.pdf",  
  plot = umap_plot,                
  width = 1150/300,                     
  height = 970/300,                      
  units = "in"                     
)
Spatial_plot<- SpatialDimPlot(seurat, group.by = "LouvainBayesian",label = TRUE,
                               label.size = 3 ,pt.size.factor = 2, alpha = 1,
                               stroke = NA)
Spatial_plot
ggsave(
    filename = "Spatial_cluster.pdf",  
    plot = Spatial_plot,                
    width = 1500/300,                     
    height = 1120/300,                      
    units = "in"                    
  )
#####-------------------part2:Comparative analysis of IBC and DCIS--------------
DCIS <- subset(seurat,idents = c(1,2,11))
IBC <- subset(seurat,idents = c(7))
IBC_vs_DCIS <- FindMarkers(seurat, ident.1 = c(7), ident.2 = c(1,2,11), min.pct = 0.25)
IBC_vs_DCIS$threshold = factor(ifelse(IBC_vs_DCIS$p_val_adj < 0.05 & abs(IBC_vs_DCIS$avg_log2FC) >= 0.5, 
                                      ifelse(IBC_vs_DCIS$avg_log2FC>= 0.5, 'Up','Down'),'NoSignifi'),
                               levels=c('Up','Down','NoSignifi'))

write.csv(IBC_vs_DCIS,file = "IBC_vs_DCIS.csv")