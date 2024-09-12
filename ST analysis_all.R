rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
devtools::load_all(".\\BayesSpace")
source("Load10X_Spatial_change.R")
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
library(ggpubr)

data.dir = ".//outs/"
seurat <- Load10X_Spatial_change(data.dir = data.dir,
                                 filename = "filtered_feature_bc_matrix",
                                 image.name = "tissue_hires_image.png")
seurat <- NormalizeData(seurat,normalization.method = "LogNormalize",scale.factor = 10000)
seurat <- FindVariableFeatures(seurat,selection.method = "vst",nfeatures = 2000)
all.genes <- VariableFeatures(seurat)
seurat <- ScaleData(seurat,features = all.genes)
seurat <- RunPCA(seurat, features = all.genes)
seurat<-JackStraw(seurat,num.replicate=10)
seurat<-ScoreJackStraw(seurat,dims=1:20)
JackStrawPlot(seurat,dims=1:20)
ndim = 15
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
set.seed(102)
sce = spatialPreprocess(sce,platform = "ST",n.PCs = ndim,n.HVGs = 2000,log.normalize = T)
q = length(unique(seurat$seurat_clusters))
set.seed(150)
sce <- spatialCluster(sce,
                      q = q, 
                      platform = "ST",
                      d = ndim, 
                      nrep = 10000, 
                      burn.in = 100,
                      init.method = "louvain",
                      model = "t",
                      gamma = 3)
seurat@meta.data = cbind(seurat@meta.data, LouvainBayes = as.factor(sce$spatial.cluster))
Idents(seurat) <- seurat@meta.data$LouvainBayes
seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:ndim)
seurat <- RunTSNE(seurat, reduction = "pca", dims = 1:ndim)
umap_plot <- DimPlot(seurat, reduction = "umap", label = TRUE)
tsne_plot <- DimPlot(seurat, reduction = "tsne", label = TRUE)
SpatialDimPlot(seurat,label.size = 3,pt.size.factor = 3,alpha = 1)

color_group <- ggplot_build(umap_plot)
colour = color_group[["data"]][[1]][["colour"]]
group =color_group[["data"]][[1]][["group"]]
info <- data.frame(cbind(colour,group))
info <- unique(arrange(info,group))
mycol <- as.character(info$colour)
mycol2 <- setNames(info$colour,as.character(info$group))




###EMT Score
gene_sets <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::filter(gs_name=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
EMT <-list(EMT_signaling=unique(gene_sets$gene_symbol)) 
seurat <- AddModuleScore(seurat,features = EMT)
colnames(seurat@meta.data)[colnames(seurat@meta.data) == "Cluster1"] <- "EMT_score"
ggplot(seurat@meta.data, aes(x=LouvainBayes, y=EMT_score,fill=LouvainBayes)) +
  geom_violin() +
  geom_boxplot(width = 0.2 ,fill = "white",outliers = FALSE)+
  scale_fill_manual(values=mycol)+
  labs(x = "", y = "EMT score") +
  theme_classic()+
  theme(legend.position = "none",
        axis.text = element_text(size=12,face="plain",color="black"))
SpatialFeaturePlot(seurat, "EMT_score",pt.size.factor = 3,alpha = 1)

###CNV analysis
coords_df <- Seurat::GetTissueCoordinates(seurat) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'barcodes')
colnames(coords_df)[2:3] <- c('y','x')
spata_obj <- SPATA2::transformSeuratToSpata(seurat_object = seurat,
                                            method = 'spatial',
                                            sample_name = 'BRCA',
                                            assay_name =  'Spatial',
                                            coords_from = coords_df)
spata_obj <- setCoordsDf(spata_obj, coords_df)
spata_obj <-
  runCnvAnalysis(
    object = spata_obj,
    directory_cnv_folder = "BRCA_cnv_obj", 
    cnv_prefix = "Chr")
cnv_results <- getCnvResults(spata_obj)
chr_cnv <- plotCnvHeatmap(object = spata_obj, across = "LouvainBayes", clrp_adjust = mycol2)
chr_cluster<- plotCnvLineplot(
  object = spata_obj,
  across = "LouvainBayes",
  n_bins_bcsp = 1000,
  n_bins_genes = 1000,
  nrow = 3
)
data_test <- as.data.frame(spata_obj@cnv[["BRCA"]][["cnv_df"]])
data_test2 <- data_test[,-1]
data_test2$mean <- rowMeans(data_test2)
barcodes <- as.data.frame(spata_obj@fdata[["BRCA"]][["barcodes"]])
names(barcodes) <- "barcodes"
data_test$mean <- data_test2$mean
result <- merge(barcodes,data_test,by = "barcodes",all.x = TRUE)
spata_obj@fdata[["BRCA"]][["mean"]] <- result$mean
cluster_bayes <- as.data.frame(spata_obj@fdata[["BRCA"]][["LouvainBayes"]])
names(cluster_bayes) <- "bayes"
data <- cbind(barcodes,cluster_bayes)
result <- merge(result,data,by = "barcodes",all.x = TRUE)
result2 <- result[,c(26,27)]
result_long <- result2 %>%
  pivot_longer(cols = c(mean), names_to = "variable", values_to = "mean")
result_long$bayes <- factor(result_long$bayes, levels = unique(result_long$bayes), ordered = FALSE) 
result_long$bayes <- factor(result_long$bayes, levels = c(1:q))
ggplot(result_long, aes(x = bayes, y = mean, fill = bayes)) +
  geom_violin() +
  geom_boxplot(width = 0.2 ,fill = "white",outliers = FALSE)+
  scale_fill_manual(values=mycol )+
  scale_y_continuous(limits = c(min(result2$mean), max(result2$mean)))+
  labs(x = "", y = "CNV result") +
  theme_classic()+
  theme(legend.position = "none",
        axis.text = element_text(size=12,face="plain",color="black"))



###### Metabolite pathway analysis
library(GSVA)
library(gson)
library(dplyr)
library(GSEABase)
library(ggplot2)
library(ggpubr)
library(reshape2)
hallmark <- read.gmt("metabolism pathway from KEGG.gmt")
hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)
matrix_all = gsva(expr = as.matrix(seurat@assays$Spatial@data), 
                  hallmark.list, 
                  kcdf="Gaussian",
                  method="ssgsea", 
                  parallel.sz=3)
matrix_all <- as.data.frame(t(matrix_all))
louvainBayes <- data.frame(seurat$LouvainBayes)
metabolite_data_result <- merge(louvainBayes,matrix_all,by = 0)
rownames(metabolite_data_result) <- metabolite_data_result[,1]
metabolite_data_result <- metabolite_data_result[,-1]
metabolite_data_result$seurat.LouvainBayes <- as.factor(metabolite_data_result$seurat.LouvainBayes)
metabolite_data_result <- metabolite_data_result %>%
  group_by(seurat.LouvainBayes) %>%
  summarise(across(where(is.numeric), ~ mean(., na.rm = TRUE)))
metabolite_data_result <- data.frame(metabolite_data_result)
rownames(metabolite_data_result) <- paste0("cluster",metabolite_data_result$seurat.LouvainBayes)
metabolite_data_result <- metabolite_data_result[,-1]
metabolite_data_result <- scale(metabolite_data_result)
metabolite_data_result <- data.frame(t(metabolite_data_result))
metabolite_data_result$Pathway <- rownames(metabolite_data_result)
rownames(metabolite_data_result) <- NULL
long_data <- melt(metabolite_data_result, varnames = c("Pathway"), value.name = "Value")
ggplot(long_data, aes(x = variable, y = Pathway, fill = Value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Heatmap of Metabolic Pathways",
       x = " ",
       y = " ",
       fill = "Value")
metabolite_data_result2 <- cbind(data,matrix_all)
metabolite_data_result2 <- metabolite_data_result2[,c("bayes","Sphingolipid metabolism")]
ggplot(metabolite_data_result2, aes(x=bayes, y=`Sphingolipid metabolism`,fill=bayes)) +
  geom_boxplot(notch=TRUE)+
  scale_fill_manual(values=mycol)+
  labs(x = "", y = "Sphingolipid metabolism") +
  theme_classic()+
  theme(legend.position = "none",
        axis.text = element_text(size=12,face="plain",color="black"))


#####Monocle Analysis
devtools::load_all(".\\monocle")
combin.data<- subset(seurat, idents = c(1,4,5,6))
expression_matrix = combin.data@assays$Spatial@counts
cell_metadata <- data.frame(group = combin.data[['orig.ident']],clusters = Idents(combin.data))
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix), stringsAsFactors = F) 
rownames(gene_annotation) <- rownames(expression_matrix)
pd <- new("AnnotatedDataFrame", data = cell_metadata)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
HSMM <- newCellDataSet(expression_matrix,
                       phenoData = pd,
                       featureData = fd,
                       expressionFamily=negbinomial.size())
HSMM_myo <- estimateSizeFactors(HSMM)
HSMM_myo <- estimateDispersions(HSMM_myo)
diff_test_res1 <- differentialGeneTest(HSMM_myo,fullModelFormulaStr = '~clusters', cores = 4)
ordering_genes <- subset(diff_test_res1, qval < 0.05)[,'gene_short_name']
HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)
HSMM_myo <- reduceDimension(HSMM_myo, max_components=2, method = 'DDRTree')
HSMM_myo <- orderCells(HSMM_myo)
plot_cell_trajectory(HSMM_myo, color_by = "clusters",cell_size =0.8)
plot_cell_trajectory(HSMM_myo,color_by = "Pseudotime",cell_size =0.8)
plot_cell_trajectory(HSMM_myo,color_by = "State",cell_size =0.8)
branch_point_choose = 2
BEAM_res <- BEAM(HSMM_myo, branch_point = branch_point_choose,progenitor_method = "duplicate",cores = 4)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
plot_genes_branched_heatmap(HSMM_myo[row.names(subset(BEAM_res,
                                                      qval < 1e-4)),],
                            branch_point = branch_point_choose,
                            num_clusters = 4,
                            cores = 4,
                            use_gene_short_name = T,
                            show_rownames = T,
                            return_heatmap = T)
top30_gene = BEAM_res[order(BEAM_res$qval),][1:30,'gene_short_name']
plot_genes_branched_heatmap(HSMM_myo[top30_gene,],
                            branch_point = branch_point_choose,
                            num_clusters = 4,
                            cores = 4,
                            use_gene_short_name = T,
                            show_rownames = T,
                            return_heatmap = T)
gene <- row.names(subset(fData(HSMM_myo),
                         gene_short_name %in% c("PSAP","CYBA")))
sample <- plot_genes_branched_pseudotime(HSMM_myo[gene,],
                               branch_point = branch_point_choose,
                               color_by = "clusters",
                               ncol = 2)
sample_data <- ggplot_build(sample)
spot_data <- sample_data$data[[1]]
names(spot_data)[1] <- "Cluster"
line_data <- sample_data$data[[2]]
spot_data[,1][spot_data[,1]=="#7CAE00"] <- "4"
spot_data[,1][spot_data[,1]=="#F8766D"] <- "1"
spot_data[,1][spot_data[,1]=="#00BFC4"] <- "5"
spot_data[,1][spot_data[,1]=="#C77CFF"] <- "6"
spot_data[,1] <- ordered(spot_data[,1], levels=c("1","4","5","6"))
spot_1 <- spot_data[spot_data$PANEL==1,]
spot_2 <- spot_data[spot_data$PANEL==2,]
line_1 <- line_data[line_data$PANEL==1,]
line_2 <- line_data[line_data$PANEL==2,]
ggplot()+
  geom_point(data = spot_1, aes(x=x, y=y, color=Cluster), size=2, alpha=0.8)+
  theme_cowplot()+
  scale_color_manual(values = c( "#F8766D","#7CAE00", "#00BFC4", "#C77CFF"))+
  geom_line(data = line_1[line_1$group==1,], aes(x=x, y=y,linetype = "Branch1"), lwd=1)+
  geom_line(data = line_1[line_1$group==2,], aes(x=x, y=y,linetype = "Branch2"), lwd=1)+ 
  labs(x = "Pseudotime(stretched)", y = "Expression")+
  scale_linetype_manual(name = " ",  
                        values = c(1, 3),  
                        labels = c("Cell fate 1", "Cell fate 2")) +
  ggtitle("CYBA") +  
  theme(plot.title = element_text(hjust = 0.5))  

ggplot()+
  geom_point(data = spot_2, aes(x=x, y=y, color=Cluster), size=2, alpha=0.8)+
  theme_cowplot()+
  scale_color_manual(values = c( "#F8766D","#00C094", "#00B6EB", "#A58AFF"))+
  geom_line(data = line_2[line_2$group==1,], aes(x=x, y=y,linetype = "Branch1"), lwd=1)+
  geom_line(data = line_2[line_2$group==2,], aes(x=x, y=y,linetype = "Branch2"), lwd=1)+ 
  labs(x = "Pseudotime(stretched)", y = "Expression")+
  scale_linetype_manual(name = " ",  
                        values = c(1, 3),  
                        labels = c("Cell fate 1", "Cell fate 2")) +
  ggtitle("PSAP") +  
  theme(plot.title = element_text(hjust = 0.5))  


#####correlation of PSAP with Sphingolipid metabolism
PSAP_gene_expression <- log2(exprs(HSMM_myo)['PSAP',]+1)
PSAP_with_metabolite_expression <- merge(PSAP_gene_expression,matrix_all,by = 0)
df <- data.frame(PSAP = PSAP_with_metabolite_expression$x,
                 Metabolite = scale(PSAP_with_metabolite_expression$`Sphingolipid metabolism`))
ggscatter(df, x = "PSAP", 
          y = "Metabolite",
          add = "reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "pearson",
          xlab = "PSAP expression",
          ylab = "Sphingolipid metabolism",
          title = "") +
  theme_classic()

save(seurat,file = "no_MIA.RData")


library(Seurat)
library(dplyr)
library(data.table)
source("Load10X_Spatial_change.R")

seurat_data <- Read10X(data.dir = "./GSE167036")
data_seurat <- CreateSeuratObject(counts = seurat_data,
                                  min.features = 200,
                                  min.cells = 3)   
scRNA <- subset(data_seurat,idents = c("CA1","CA2","CA3","CA4","CA5","CA6","CA7","CA8"))
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 3000)
all.genes <- VariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, features = all.genes)
scRNA <- RunPCA(scRNA, features = all.genes)
cell_type <- read.csv(file = "GSE167036_meta_all_change.csv",row.names = 1)
cell_type_select <- cell_type[cell_type$orig.ident %in% c("CA1", "CA2","CA3","CA4","CA5","CA6","CA7","CA8"), ]
scRNA@meta.data = cbind(scRNA@meta.data, cell_type = as.factor(cell_type_select$main_label)) 
Idents(scRNA) <- scRNA@meta.data$cell_type
set.seed(123)
maintype_marker=FindAllMarkers(scRNA,logfc.threshold = 0.5,only.pos = T)
maintype_marker=maintype_marker%>%filter(p_val_adj < 1e-05)
maintype_marker$d=maintype_marker$pct.1 - maintype_marker$pct.2
maintype_marker=maintype_marker%>%filter(d > 0.2)
maintype_marker=maintype_marker%>%arrange(cluster,desc(avg_log2FC))
maintype_marker=as.data.frame(maintype_marker)
save(scRNA,maintype_marker,file = "MIA-scRNA.RData")


library(Seurat)
library(dplyr)
library(data.table)
library(ggplot2)
library(spacexr)
library(STdeconvolve)
library(ggsci)
library(RColorBrewer)
library(tidyr)

load("MIA-scRNA.RData")
load("no_MIA.RData")

####MIA Analysis
region_marker=FindAllMarkers(seurat,logfc.threshold = 0,only.pos = T)
region_marker=region_marker%>%filter(p_val_adj < 0.1)
region_marker$d=region_marker$pct.1 - region_marker$pct.2
region_marker=region_marker%>%filter(d > 0.05)
region_marker=region_marker%>%arrange(cluster,desc(avg_log2FC))
region_marker=as.data.frame(region_marker)
st.clusts <- Idents(seurat) %>% levels()
N <- length(st.clusts)
st.marker.list <- vector(mode = 'list', length = N)
names(st.marker.list) <- st.clusts
for(i in st.clusts) {
  st.marker.list[[i]] <- region_marker[region_marker$cluster == i,'gene']
}
sc.clusts <- Idents(scRNA) %>% levels()
M <- length(sc.clusts)
sc.marker.list <- vector(mode = 'list', length = M)
names(sc.marker.list) <- sc.clusts
for (i in sc.clusts) {
  sc.marker.list[[i]] <- maintype_marker[maintype_marker$cluster == i,'gene']
}
N <- length(st.clusts) ; M <- length(sc.clusts)
MIA.results <- matrix(0,nrow = M, ncol = N)
row.names(MIA.results) <- sc.clusts
colnames(MIA.results) <- st.clusts
gene.universe <- length(rownames(seurat))
for (i in 1:N) {
  for (j in 1:M) {
    genes1 <- st.marker.list[[st.clusts[i]]]
    genes2 <- sc.marker.list[[sc.clusts[j]]]
    A <- length(intersect(genes1,genes2))
    B <- length(genes1)
    C <- length(genes2)
    enr <- -log10(phyper(A, B, gene.universe-B, C, lower.tail = FALSE))
    dep <- -log10(1-phyper(A, B, gene.universe-B, C, lower.tail = FALSE))
    if (enr < dep) {
      MIA.results[j,i] = -dep
    } else {
      MIA.results[j,i] = enr
    }
  }
}
MIA.results[is.infinite(MIA.results)] <- 0
heatmap_df <- data.frame('Cell_types' = reshape2::melt(MIA.results)[,1],
                         'Cluster' = reshape2::melt(MIA.results)[,2],
                         enrichment = reshape2::melt(MIA.results)[,3])
heatmap_df$Cluster <- as.factor(heatmap_df$Cluster)
ggplot(data = heatmap_df, aes(x = Cluster, y = Cell_types, fill = enrichment)) +
  geom_tile() +
  scale_fill_gradient2(low = "navyblue", high = "red", mid = "white",
                       midpoint = 0, limit = c(0,110), space = "Lab", oob=scales::squish,
                       name="Enrichment \n -log10(p)") +
  scale_x_discrete(labels = function(x) as.character(x))+
  ylim(heatmap_df$Cell_types %>% levels() %>% sort() %>% rev()) +
  theme_minimal()


####RCTD Analysis
Idents(scRNA) <- "celltype"
counts <- scRNA[["RNA"]]$counts
cluster <- as.factor(scRNA$cell_type)
names(cluster) <- colnames(scRNA)
nUMI <- scRNA$nCount_RNA
names(nUMI) <- colnames(scRNA)
set.seed(123)
reference <- Reference(counts, cluster, nUMI)

counts <- seurat[["Spatial"]]$counts
coords <- GetTissueCoordinates(seurat)
colnames(coords) <- c("x", "y")
coords[is.na(colnames(coords))] <- NULL
query <- SpatialRNA(coords, counts, colSums(counts))


RCTD <- create.RCTD(query, reference, max_cores = 1)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")
barcodes <- colnames(RCTD@spatialRNA@counts)
weights <- RCTD@results$weights
norm_weights <- normalize_weights(weights)
plt <- vizAllTopics(theta = as.matrix(norm_weights),
                    pos = coords,
                    topicOrder=seq(ncol(norm_weights)),
                    topicCols=rainbow(ncol(norm_weights)),
                    groups = NA,
                    group_cols = NA,
                    r = 20, # size of scatterpies; adjust depending on the coordinates of the pixels
                    lwd = 0.3,
                    showLegend = TRUE,
                    plotTitle = "")
plt <- plt + ggplot2::guides(fill=ggplot2::guide_legend(ncol=1))
plt



expression_matrix <- data.frame(as.matrix(norm_weights))
cluster_result <- data.frame(seurat$LouvainBayes)
result <- merge(cluster_result,expression_matrix,by = 0)
rownames(result) <- result$Row.names
result <- result[,-1]
result$seurat.LouvainBayes <- as.factor(result$seurat.LouvainBayes)
mean_values <- result %>%
  group_by(seurat.LouvainBayes) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))
mean_values <- mean_values %>%
  pivot_longer(
    cols = -seurat.LouvainBayes,  
    names_to = "cell_type",  
    values_to = "percent"  
  )


p1<- ggplot(data=mean_values,aes(x=seurat.LouvainBayes,y=percent,fill=cell_type))+
  geom_bar(stat = "identity", 
           position="fill", 
           width = 0.8,   
           color="black", 
           linewidth=0.3)+ 
  theme_bw()+
  scale_y_continuous(expand = expansion(mult=c(0,0.05),add=c(0,0)))+   
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5))

p1





