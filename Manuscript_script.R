##load required packages
ipak <- function(pkg){new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if(length(new.pkg))
  install.packages(new.pkg, dependencies = TRUE)
sapply(pkg, require, character.only = TRUE)}
packages <- c("Seurat", "dplyr", "ggplot2", "tidyverse", "patchwork", "devtools", "stringi",
              "reticulate", "rlang", "caTools", "cowplot", "MAST", "viridis",
              "Matrix", "RCurl","SingleCellExperiment")
ipak(packages)

#Seurat pipeline was mainly used for data analysis
#1.load 5' 10x single cell sequence data (P can be replaced with P1~P11 and N1~N3),create count matrix and add metadata.
#2.create metadata(nCount_RNA, nFeature_RNA, percent.mt, mitoratio, log10GenePerUMI, CDR3sequence, clonotype_id) 
P_5P.data <- Read10X_h5("~/project/2021_CTCL_data/5P_Tcell/P-T_filtered_feature_bc_matrix.h5")
P_5P <-CreateSeuratObject(counts = P_5P.data, project = "P_5P_tcr", min.features = 100)
##add paired TCR sequencing data as seurat object metadata
P_tcr <- read.csv("~/project/2021_CTCL_data/5P_TCR/Filtered_contig_annotations/P_filtered_contig_annotations.csv", sep = ",")
P_tcr <- P_tcr[!duplicated(P_tcr$barcode), ]
P_tcr <- P_tcr[, c("barcode", "raw_clonotype_id")]
names(P_tcr)[names(P_tcr) == "raw_clonotype_id"] <- "clonotype_id"
P_clono <- read.csv("~/project/2021_CTCL_data/5P_TCR/clonotypes/P_clonotypes.csv", sep = ",")
P_tcr <- merge(P_tcr, P_clono[, c("clonotype_id", "cdr3s_aa")])
P_tcr <- P_tcr[, c(2,1,3)]
rownames(P_tcr) <- P_tcr[, 1]
P_tcr[, 1] <- NULL
P_merged_T <- AddMetaData(object = P_5P, metadata = P_tcr)
##add percent.mt, mitoRatio,log10GenePerUMI into metadata
P_merged_T$percent.mt <- PercentageFeatureSet(P_merged_T, pattern = "^MT-")
P_merged_T$mitoRatio <- P_merged_T@meta.data$percent.mt / 100
P_merged_T$log10GenePerUMI <- log10(P_merged_T$nFeature_RNA) / log10(P_merged_T$nCount_RNA)
metadata <- P_merged_T@meta.data

saveRDS(P_merged_T, file = file.path("~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/data/", "P_merged_T.rds"))
P_merged_T <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/data/P_merged_T.rds")

#filter low quality cells and genes
# QC visualization
P_QC_plot1 <-FeatureScatter(P_merged_T, feature1 = "nFeature_RNA", feature2 = "percent.mt", pt.size = 0.5)
P_QC_plot1

#filter low quality cells
filtered_P_merged_T <- subset(x=P_merged_T,
                                   subset= (nCount_RNA >=500) &
                                     (nFeature_RNA >=250) &
                                     (log10GenePerUMI>0.8) &
                                     (mitoRatio <0.2))

Filtered_P_merged_T  ###(33538 features across 12568 samples within 1 assay )

#filter genes
##remove TCR genes 
counts <- GetAssayData(filtered_P_merged_T, assay = "RNA")
remain <- setdiff(1:nrow(counts), grep("^TRAV|^TRAJ|^TRBV|^TRBJ|^TRBD|^TRDD|^TRDJ|^TRDV|^TRGV|^TRGJ", rownames(counts)))
remianed_counts <- counts[remain, ]
filtered_P_merged_T <- subset(filtered_JA1819_merged_T, features = rownames(remianed_counts))
filtered_P_merged_T  ##33311 features across 12568 samples within 1 assay  

##remove "0" gene
counts <- GetAssayData(object = filtered_P_merged_T, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix:: rowSums(nonzero) >=10
filtered_counts <- counts[keep_genes, ]
filtered_P_merged_T <- CreateSeuratObject(filtered_counts, meta.data = filtered_JA1819_merged_T@meta.data)

filtered_P_merged_T ##(14960 features across 12568 samples within 1 assay )


#recheck the QC metrics
filtered_P_merged_T@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) 
filtered_P_merged_T_plot <-FeatureScatter(filtered_P_merged_T, feature1 = "nFeature_RNA", feature2 = "percent.mt", pt.size = 0.5)
filtered_P_merged_T_plot
##save data
saveRDS(filtered_P_merged_T, file = file.path("~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/data/", "filtered_P_merged_T.rds"))
filtered_P_merged_T <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/data/filtered_P_merged_T.rds")


#filter non CD4+ T cells via removing non-CD4 T cell clusters
filtered_P_merged_T <-NormalizeData(filtered_P_merged_T, verbose = TRUE) 
filtered_P_merged_T <- FindVariableFeatures(filtered_P_merged_T, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(filtered_P_merged_T)
filtered_P_merged_T <- ScaleData(filtered_P_merged_T, features = all.genes)
filtered_P_merged_T <- RunPCA(filtered_P_merged_T, features = VariableFeatures(object = filtered_P_merged_T), verbose = TRUE)
DimPlot(filtered_P_merged_T, reduction = "pca")
DimHeatmap(filtered_P_merged_T, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(filtered_P_merged_T)
filtered_P_merged_T <- FindNeighbors(filtered_P_merged_T, dims = 1:20, verbose = TRUE)
filtered_P_merged_T <- FindClusters(filtered_P_merged_T, resolution = 0.02, verbose = TRUE)
filtered_P_merged_T <- RunUMAP(filtered_P_merged_T, dims = 1: 20)
DimPlot(filtered_P_merged_T, reduction = "umap", label = TRUE, label.size = 6, pt.size = 0.8)

#Use dominant clonotype to annotate the CTCL cluster
clonotype1 <- subset(filtered_P_merged_T, subset = clonotype_id == "clonotype1")
clonotype1 <- Cells(clonotype1)
clonotype_dominant_plot <- DimPlot(filtered_P_merged_T, cells.highlight = clonotype1, reduction = "umap") 
clonotype_dominant_plot

#recheck the quality of the cells in each cluster
VlnPlot(filtered_P_merged_T, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

#Annotate populations of each cluster
FeaturePlot(filtered_P_merged_T, features = c("CD3E", "CD3D","CD8A","CD8B", "CD4",
                                                   "ITGB1", "CD79A", "NKG7","GNLY", "MS4A1", 
                                                   "CD79A", "MKI67", "TOP2A" ))

#subset normal CD4 T cell cluster
n_CD4_P <- subset(filtered_P_merged_T, ident= c(1, 3, 4))
n_CD4_P <-NormalizeData(n_CD4_P, verbose = TRUE) 
n_CD4_P <- FindVariableFeatures(n_CD4_P, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(n_CD4_P)
n_CD4_P <- ScaleData(n_CD4_P, features = all.genes)
n_CD4_P <- RunPCA(n_CD4_P, features = VariableFeatures(object = n_CD4_P), verbose = TRUE)
DimPlot(n_CD4_P, reduction = "pca")
DimHeatmap(n_CD4_P, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(n_CD4_P)
n_CD4_P <- FindNeighbors(n_CD4_P, dims = 1:15, verbose = TRUE)
n_CD4_P <- FindClusters(n_CD4_P, resolution = 0.5, verbose = TRUE)
n_CD4_P <- RunUMAP(n_CD4_P, dims = 1: 15)
DimPlot(n_CD4_P, reduction = "umap", label = TRUE, label.size = 6, pt.size = 0.8)
#Find tumor clonotype
clonotype1 <- subset(n_CD4_P, subset = clonotype_id == "clonotype1")
clonotype1 <- Cells(clonotype1)
clonotype_dominant_plot <- DimPlot(n_CD4_P, cells.highlight = clonotype1, reduction = "umap") 
clonotype_dominant_plot
VlnPlot(n_CD4_P, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
FeaturePlot(n_CD4_P, features = c("CD3E", "CD3D","CD8A","CD8B", "CD4",
                                       "ITGB1","CD14", "CD79A", "NKG7","GNLY", "MS4A1", 
                                       "CD79A", "MKI67", "TOP2A", "PCLAF","TYMS", "TRGC", "TRDC", "TRAC", "TRBC" ))
n_CD4_P <- subset(n_CD4_P, ident = c(0,1,2,3,4,6))
saveRDS(n_CD4_P, file = file.path("~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/data/", "n_CD4_JA1819.rds"))

#subset CTCL cluster
CTCL_P <- subset(filtered_P_merged_T, ident= c(0,2))
CTCL_P <-NormalizeData(CTCL_P, verbose = TRUE) 
CTCL_P <- FindVariableFeatures(CTCL_P, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(CTCL_P)
CTCL_P <- ScaleData(CTCL_P, features = all.genes)
CTCL_P <- RunPCA(CTCL_P, features = VariableFeatures(object = CTCL_P), verbose = TRUE)
DimPlot(CTCL_P, reduction = "pca")
DimHeatmap(CTCL_P, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(CTCL_P)
CTCL_P <- FindNeighbors(CTCL_P, dims = 1:15, verbose = TRUE)
CTCL_P <- FindClusters(CTCL_P, resolution = 0.5, verbose = TRUE)
CTCL_P <- RunUMAP(CTCL_P, dims = 1: 15)
DimPlot(CTCL_P, reduction = "umap", label = TRUE, label.size = 6, pt.size = 0.8)
#Find tumor clonotype
clonotype1 <- subset(CTCL_P, subset = clonotype_id == "clonotype1")
clonotype1 <- Cells(clonotype1)
clonotype_dominant_plot <- DimPlot(CTCL_P, cells.highlight = clonotype1, reduction = "umap") 
clonotype_dominant_plot
VlnPlot(CTCL_P, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
FeaturePlot(CTCL_P, features = c("CD3E", "CD3D","CD8A","CD8B", "CD4",
                                      "ITGB1","NKG7","GNLY", "MS4A1", 
                                      "CD79A", "MKI67", "TOP2A", "PCLAF","TYMS" ))

CTCL_P <- subset(CTCL_P, ident= c(0,1,2,3,4,5))
saveRDS(CTCL_P, file = file.path("~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/data/", "CTCL_P.rds"))


#final qualified cells and genes
n_CD4_P <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/data/n_CD4_P.rds")
CTCL_P <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/data/CTCL_P.rds")
P_5P_finalQC <- merge(x= n_CD4_P, y= CTCL_P)
P_5P_finalQC
P_5P_finalQC$RNA_snn_res.0.02 <- NULL
P_5P_finalQC$seurat_clusters <- NULL
P_5P_finalQC$RNA_snn_res.0.5 <- NULL

saveRDS(P_5P_finalQC, file = file.path("~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/data/", "P_5P_finalQC.rds"))
P_5P_finalQC <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/data/P_5P_finalQC.rds" )




##Single sample analysis
P_5P_U <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/data/P_5P_finalQC.rds")

P_5P_U <-NormalizeData(P_5P_U, verbose = TRUE) 
P_5P_U <- FindVariableFeatures(P_5P_U, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(P_5P_U)
P_5P_U <- ScaleData(P_5P_U, features = all.genes)
P_5P_U <- RunPCA(P_5P_U, features = VariableFeatures(object = P_5P_U), verbose = TRUE)
DimPlot(P_5P_U, reduction = "pca")
DimHeatmap(P_5P_U, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(P_5P_U)
P_5P_U <- FindNeighbors(P_5P_U, dims = 1:20, verbose = TRUE)
P_5P_U <- FindClusters(P_5P_U, resolution = 0.01, verbose = TRUE)
P_5P_U <- RunTSNE(P_5P_U, dims = 1: 20)
DimPlot(P_5P_U, reduction = "tsne", label = TRUE, label.size = 6, pt.size = 0.8)

# label cell type of clusters
new.cluster.ids <- c("CTCL", "Normal_CD4T")
names(new.cluster.ids) <- levels(P_5P_U)
P_5P_U <- RenameIdents(JA1819_5P_U, new.cluster.ids)
DimPlot(P_5P_U, label = TRUE, label.size = 8, pt.size = 0.5, reduction = "tsne") + NoLegend()
#subset normal CD4 T 
normal_P <- subset(P_5P_U, ident = "Normal_CD4T")
normal_P <-NormalizeData(normal_P, verbose = TRUE) 
normal_P <- FindVariableFeatures(normal_P, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(normal_P)
normal_P <- ScaleData(normal_P, features = all.genes)
normal_P <- RunPCA(normal_P, features = VariableFeatures(object = normal_P), verbose = TRUE)
DimPlot(normal_P, reduction = "pca")
ElbowPlot(normal_P)
normal_P <- FindNeighbors(normal_P, dims = 1:15, verbose = TRUE)
normal_P <- FindClusters(normal_P, resolution = 0.1, verbose = TRUE)
normal_P <- RunTSNE(normal_P, dims = 1: 15)
DimPlot(normal_P, label = TRUE, label.size = 6, pt.size = 0.8,reduction = "tsne")
clonotype1 <- subset(normal_P, subset = clonotype_id == "clonotype1")
clonotype1 <- Cells(clonotype1)
clonotype_dominant_plot <- DimPlot(normal_P, cells.highlight = clonotype1, reduction = "tsne") 
clonotype_dominant_plot

new.cluster.ids <- c("Non-CTCL CD4+ T cell population 1","Non-CTCL CD4+ T cell population 2")
names(new.cluster.ids) <- levels(normal_P)
normal_P <- RenameIdents(normal_P, new.cluster.ids)
DimPlot(normal_P, label = TRUE, label.size = 8, pt.size = 0.5,reduction = "tsne" ) + NoLegend()
saveRDS(normal_P, file = file.path("~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Phate_subpopulations/rsd_file/", "normal_P.rds"))
normal_P <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Phate_subpopulations/rsd_file/normal_P.rds")

#subset tumor cluster
tumor_P <- subset(P_5P_U, ident = "CTCL")
## integrate subseted idents with original UMAP object 
P_5P_U$subseted_cluster <- as.character(Idents(P_5P_U))
P_5P_U$subseted_cluster[Cells(normal_P)] <- paste("Normal_CD4T", Idents(normal_P))
P_5P_U$subseted_cluster[Cells(tumor_P)] <- paste("CTCL", Idents(tumor_P))
DimPlot(P_5P_U, group.by = "subseted_cluster",reduction = "tsne" )


saveRDS(P_5P_U, file = file.path("~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Phate_subpopulations/rsd_file/", "P_phate.rds"))
P_phate <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Phate_subpopulations/rsd_file/P_phate.rds" )

##PHATE trajectory analysis Figure2
#load/install required packages
devtools::install_github("scottgigante/seurat", ref="patch/add-PHATE-again")  # Installation of Seurat with PHATE function embedded
library(Seurat)
P_phate <- readRDS("~/Downloads/P_phate.rds")
P_phate <- RunPHATE(P_phate, knn = 5L, decay = 40L)
Data.S <- P_phate
Data.S <- RunUMAP(P_phate, dims = 1: 20)
clonotype1 <- subset(Data.S, subset = clonotype_id == "clonotype1")  ##dominant clonotype plotting on the phate
clonotype1 <- Cells(clonotype1)
DimPlot(Data.S, reduction = "umap",label.size = 6, pt.size = 1, cells.highlight = clonotype1)+ scale_color_manual(labels = c("Other_clone", "CTCL_clone"), values = c("#808080", "#FF8C00"))+ NoLegend()
DimPlot(object = Data.S,group.by = "subseted_cluster", label = FALSE, reduction = 'tsne', pt.size = 0.5,cols = c("CTCL CTCL" = "#FF0000", "Normal_CD4T normal_clonotype1-"= "#0000FF", "Normal_CD4T normal_clonotype1+" = "#8B008B"))+NoLegend() #tsne plot
DimPlot(Data.S, label = FALSE, reduction = "phate", group.by = "subseted_cluster",pt.size = 0.5,cols = c("CTCL CTCL" = "#FF0000", "Normal_CD4T normal_clonotype1-"= "#0000FF", "Normal_CD4T normal_clonotype1+" = "#8B008B"))+NoLegend()  #original phate plot
clonotype1 <- subset(Data.S, subset = clonotype_id == "clonotype1")  ##dominant clonotype plotting on the phate
clonotype1 <- Cells(clonotype1)
DimPlot(Data.S, reduction = "phate",label.size = 6, pt.size = 1, cells.highlight = clonotype1)+ scale_color_manual(labels = c("Other_clone", "CTCL_clone"), values = c("#808080", "#FF8C00"))+ NoLegend()


##Integration analysis of 11 patients and 3 healthy controls
##Integration of 3 healthy controls 
#loading filtered seurat objects
ss_NS_merged_modified <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/data/NS_5P_finalQC.rds")
ss_N1_merged_modified <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/data/N1_5P_finalQC.rds")
ss_N2_merged_modified <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/data/N2_5P_finalQC.rds")

#load the data and add metadata
ss_NS_merged_modified@meta.data[,"category"] <- "NS"
ss_N1_merged_modified@meta.data[,"category"] <- "N1"
ss_N2_merged_modified@meta.data[,"category"] <- "N2"

#Integration
sample_combined = merge(x = ss_NS_merged_modified, y = c( ss_N1_merged_modified, ss_N2_merged_modified), add.cell.ids = c("NS", "N1", "N2"), Project = "catogory")
sample.list <- SplitObject(sample_combined, split.by = "category")
reference.list <- sample.list[c("NS", "N1", "N2")]
sample.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:15)
saveRDS(sample.anchors, file = file.path("~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/UMAP/", "sample.anchors_3_normals.rds"))
sample.anchors <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/UMAP/sample.anchors_3_normals.rds")

#UMAP and visualization
##UMAP
sample.integrated.3_normals <- IntegrateData(anchorset = sample.anchors, dims = 1:15)
DefaultAssay(sample.integrated.3_normals) <- "integrated"
sample.integrated.3_normals <- ScaleData(sample.integrated.3_normals, verbose = TRUE)
sample.integrated.3_normals <- RunPCA(sample.integrated.3_normals, npcs = 30, verbose = TRUE)
sample.integrated.3_normals <- RunUMAP(sample.integrated.3_normals, reduction = "pca", dims = 1:15)
sample.integrated.3_normals <- FindNeighbors(sample.integrated.3_normals, reduction = "pca", dims = 1:15)
sample.integrated.3_normals <- FindClusters(sample.integrated.3_normals, resolution = 0.1)
P1 <- DimPlot(sample.integrated.3_normals, reduction = "umap", group.by = "category")
plot(P1)
P2 <- DimPlot(sample.integrated.3_normals, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
plot(P2)

##visualize side-by-side##
DimPlot(sample.integrated.3_normals, reduction = "umap", split.by = "category", label = TRUE)



#Integration of 11 patients
##load the 10x data and creat seurat objects
ss_P11_merged_modified <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/UMAP/P11_5P_U_All.rds")
ss_P10_merged_modified <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/UMAP/P10_5P_U_All.rds")
ss_P9_merged_modified <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/UMAP/P9_5P_U_All.rds")
ss_P8_merged_modified <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/UMAP/P8_5P_U_All.rds")
ss_P7_merged_modified <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/UMAP/P7_5P_U_All.rds")
ss_P6_merged_modified <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/UMAP/P6_5P_U_All.rds")
ss_P5_merged_modified <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/UMAP/P5_5P_U_All.rds")
ss_P4_merged_modified <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/UMAP/P4_5P_U_All.rds")
ss_P3_merged_modified <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/UMAP/P3_5P_U_All.rds")
ss_P2_merged_modified <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/UMAP/P2_5P_U_All.rds")
ss_P1_merged_modified <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/UMAP/P1_5P_U_All.rds")

##load the data and add metadata

ss_P11_merged_modified@meta.data[,"category"] <- "P11"
ss_P10_merged_modified@meta.data[,"category"] <- "P10"
ss_P9_merged_modified@meta.data[,"category"] <- "P9"
ss_P8_merged_modified@meta.data[,"category"] <- "P8"
ss_P7_merged_modified@meta.data[,"category"] <- "P7"
ss_P6_merged_modified@meta.data[,"category"] <- "P6"
ss_P5_merged_modified@meta.data[,"category"] <- "P5"
ss_P4_merged_modified@meta.data[,"category"] <- "P4"
ss_P3_merged_modified@meta.data[,"category"] <- "P3"
ss_P2_merged_modified@meta.data[,"category"] <- "P2"
ss_P1_merged_modified@meta.data[,"category"] <- "P1"

##Integration

sample_combined = merge(x = ss_SEOG_merged_modified, y = c(ss_A14_merged_modified, ss_A23_merged_modified, ss_AU31_merged_modified, ss_AUG1914_merged_modified, ss_FEB112020_2_merged_modified, ss_FEB192020_merged_modified, ss_JA1819_merged_modified, ss_JUL18_merged_modified, ss_MA1319_merged_modified, ss_NOV1819_merged_modified), add.cell.ids = c("SEOG", "A14", "A23", "AU31", "AUG1914", "FEB112020_2", "FEB192020", "JA1819", "JUL18", "MA1319", "NOV1819"), Project = "catogory")
sample.list <- SplitObject(sample_combined, split.by = "category")
reference.list <- sample.list[c("P11", "P10", "P9", "P8", "P7", "P6", "P5", "P4", "P3", "P2", "P1")]
sample.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:15)
saveRDS(sample.anchors, file = file.path("~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/UMAP", "sample.anchors_11_patients.rds"))
sample.anchors <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/UMAP/sample.anchors_11_patients.rds")

##UMAP and visualization

##UMAP
sample.integrated.11.patients <- IntegrateData(anchorset = sample.anchors, dims = 1:15)
DefaultAssay(sample.integrated.11.patients) <- "integrated"
sample.integrated.11.patients <- FindVariableFeatures(sample.integrated.11.patients)
sample.integrated.11.patients <- ScaleData(sample.integrated.11.patients, verbose = TRUE)
sample.integrated.11.patients <- RunPCA(sample.integrated.11.patients, npcs = 30, verbose = TRUE)
sample.integrated.11.patients <- RunUMAP(sample.integrated.11.patients, reduction = "pca", dims = 1:15)
sample.integrated.11.patients <- FindNeighbors(sample.integrated.11.patients, reduction = "pca", dims = 1:15)
sample.integrated.11.patients <- FindClusters(sample.integrated.11.patients, resolution = 0.05)
P1 <- DimPlot(sample.integrated.11.patients, reduction = "umap", group.by = "category")
plot(P1)
P2 <- DimPlot(sample.integrated.11.patients, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
plot(P2)
plot_grid(P1, P2)


##visualize side-by-side##
DimPlot(sample.integrated.11.patients, reduction = "umap", split.by = "category", label = TRUE)

##Assign clonotype back to cluster and indentify tumor clusters
clonotype1 <- subset(sample.integrated.11.patients, subset = clonotype_id == "clonotype1")
clonotype1 <- Cells(clonotype1)
clonotype_dominant_plot <- DimPlot(sample.integrated.11.patients, cells.highlight = clonotype1, reduction = "umap", split.by = "category", label = TRUE,pt.size = 0.02 )
clonotype_dominant_plot

#Integration of Healhty controls and Patients together

sample.integrated.11.patients <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/UMAP/sample.integrated.11.patients.rds")
sample.integrated.3_normals <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/UMAP/sample.integrated.3_normals.rds")
DefaultAssay(sample.integrated.11.patients) <- "RNA"
DefaultAssay(sample.integrated.3_normals) <- "RNA"
sample.integrated.11.patients@meta.data[, "category"] <- "CTCL_Patients"
sample.integrated.3_normals@meta.data[, "category"] <- "Healthy_Control"

sample.list = merge(x= sample.integrated.11.patients, y = sample.integrated.3_normals, add.cell.ids = c("CTCL_Patients", "Healthy_Control"), project = "category")
sample.list <- SplitObject(sample.list, split.by = "category")
reference.list <- sample.list[c("CTCL_Patients", "Healthy_Control")]
sample.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:15)
saveRDS(sample.anchors, file = file.path("~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/UMAP/", "sample.anchors_11and 3_patients.rds"))
sample.anchors <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/UMAP/sample.anchors_11and 3_patients.rds")

##UMAP and visualization

sample.integrated.11_3 <- IntegrateData(anchorset = sample.anchors, dims = 1:15)
DefaultAssay(sample.integrated.11_3) <- "integrated"
sample.integrated.11.patients <- FindVariableFeatures(sample.integrated.11_3)
sample.integrated.11_3 <- ScaleData(sample.integrated.11_3, verbose = TRUE)
sample.integrated.11_3 <- RunPCA(sample.integrated.11_3, npcs = 30, verbose = TRUE)
sample.integrated.11_3 <- RunUMAP(sample.integrated.11_3, reduction = "pca", dims = 1:15)
sample.integrated.11_3 <- FindNeighbors(sample.integrated.11_3, reduction = "pca", dims = 1:15)
sample.integrated.11_3 <- FindClusters(sample.integrated.11_3, resolution = 0.05)
P1 <- DimPlot(sample.integrated.11_3, reduction = "umap", group.by = "category")
plot(P1)
P2 <- DimPlot(sample.integrated.11_3, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
plot(P2)




###############################
#Figure 1B
library(ggplot2)
library(scales)

clonotype_table <- as.data.frame(
  read_csv("/data/project_202108/metadata_sample.integrated.11_3.csv")
)
colnames(clonotype_table)
rownames(clonotype_table) <- clonotype_table[, 1]
clonotype_table <- clonotype_table[, -1]
head(clonotype_table)
samples <- unique(clonotype_table$orig.ident)

#1. count % of NA, single, double, triple and fouth chains of clonotype in each sample
clonotype_table$clonotype_category <- ""
clonotype_table$clonotype_category[which(is.na(clonotype_table$clonotype_id))] <- "NA"
clonotype_table$cdr3s_aa[1:10]
clonotype_table$N_chain <- stringr::str_count(clonotype_table$cdr3s_aa, ";") + 1
clonotype_table$clonotype_category[which(clonotype_table$N_chain == 1)] <- "single chain"
clonotype_table$clonotype_category[which(clonotype_table$N_chain == 2)] <- "double chain"
clonotype_table$clonotype_category[which(clonotype_table$N_chain == 3)] <- "triple chain"
clonotype_table$clonotype_category[which(clonotype_table$N_chain == 4)] <- "quadruple chain"
summary_table <- table(clonotype_table$orig.ident, clonotype_table$clonotype_category)
summary_table_normalized <- summary_table/as.numeric(table(clonotype_table$orig.ident))
openxlsx::write.xlsx(summary_table, "summary_table.xlsx")
openxlsx::write.xlsx(summary_table_normalized, "summary_table_normalized.xlsx")

#2. remove NA from clonotype_id
clonotype_table_all <- clonotype_table
clonotype_table <- clonotype_table[which(!is.na(clonotype_table$clonotype_id)), ]

#3. calculate the % of each clonotype in each sample and save as an excel file for future tacking, including ident, clonotype_id, freq%, and CDR3_aa sequence.
for (i in samples){
  print(i)
  sub_clonotype_table <- clonotype_table[which(clonotype_table$orig.ident == i), ]
  clonotype <- paste0("clonotype", 1:length(unique(sub_clonotype_table$clonotype_id)))
  freq <- as.numeric(table(sub_clonotype_table$clonotype_id)[clonotype]/nrow(sub_clonotype_table))
  CDR3_aa <- rep("", length(clonotype))
  for (j in 1:length(CDR3_aa)){
    CDR3_aa[j] <- paste(unique(sub_clonotype_table[which(sub_clonotype_table$clonotype_id == clonotype[j]), ]$cdr3s_aa), collapse = " | ")
  }
  df <- data.frame(clonotype_id = clonotype, freq = freq, cdr3s_aa = CDR3_aa)
  openxlsx::write.xlsx(df, paste0(i, "_summary.xlsx"))
}


unique(clonotype_table[which(clonotype_table$clonotype_id == "clonotype1"), ]$cdr3s_aa)
clonotype_table$clonotype1 <- "other"
for (sample_i in samples){
  sub_clonotype_table <- clonotype_table[which(clonotype_table$orig.ident == sample_i), ]
  clonotype1_cdr3s_aa <- unique(sub_clonotype_table[which(sub_clonotype_table$clonotype_id == "clonotype1"), ]$cdr3s_aa)
  print(clonotype1_cdr3s_aa)
  chains <- strsplit(clonotype1_cdr3s_aa, ";")[[1]]
  TRA_aa <- chains[grep("TRA", chains)]
  TRB_aa <- chains[grep("TRB", chains)]
  for (i in 1:length(sub_clonotype_table$clonotype1)){
    clonotype1_cdr3s_aa <- sub_clonotype_table$cdr3s_aa[i]
    chains <- strsplit(clonotype1_cdr3s_aa, ";")[[1]]
    tra_aa <- chains[grep("TRA", chains)]
    trb_aa <- chains[grep("TRB", chains)]
    if (all(tra_aa %in% TRA_aa) | all(trb_aa %in% TRB_aa)){
      sub_clonotype_table$clonotype1[i] <- "clonotype1-like"
    }
  }
  sub_clonotype_table$clonotype1[which(sub_clonotype_table$clonotype_id == "clonotype1")] <- "clonotype1"
  if (sample_i == "FEB112020_2_5P_tcr") {
    sub_clonotype_table$clonotype1[which(sub_clonotype_table$clonotype_id == "clonotype2")] <- "clonotype1"
  }
  clonotype_table[rownames(sub_clonotype_table),]$clonotype1 <- sub_clonotype_table$clonotype1
}
clonotype_table$clonotype1 <- factor(clonotype_table$clonotype1,
                                     levels = c("clonotype1", "clonotype1-like", "other"))
clonotype_table$orig.ident <- factor(clonotype_table$orig.ident,
                                     levels = c(samples[1:11], samples[12:14]))
table(clonotype_table$orig.ident, clonotype_table$clonotype1)

p <- ggplot(clonotype_table, aes(x=orig.ident, fill=clonotype1)) + geom_bar(position = "fill") + #position = "fill"
  theme(text = element_text(size=15))  + 
  scale_fill_manual(values = c("darkred", "pink", "lightgray")) + theme(legend.title = element_blank()) +
  coord_flip() + 
  scale_y_reverse(breaks=c(0,0.25,0.5,0.75,1),
                  labels = rev(c("0%","25%","50%","75%","100%"))) + 
  labs(x = "", y = "Percentage") #+ RotatedAxis() 
p



#subset CTCL_patient cluster
sample.integrated.11_3 <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/UMAP/sample.integrated.11_3.rds")
gene_integrated <- sample.integrated.11_3
gene_integrated$new_ident <- paste(gene_integrated$category, gene_integrated$seurat_clusters, sep="_")
gene_integrated$phate_cluster <- "Healthy_Normal"
gene_integrated$phate_cluster[which(gene_integrated$new_ident %in% "CTCL_Patients_0")] <- "CTCL_Normal"
gene_integrated$phate_cluster[which(gene_integrated$new_ident %in% "CTCL_Patients_12")] <- "CTCL_Normal"
gene_integrated$phate_cluster[which(gene_integrated$new_ident %in% "CTCL_Patients_13")] <- "CTCL_Normal"
gene_integrated$phate_cluster[which(gene_integrated$new_ident %in% "CTCL_Patients_1")] <- "CTCL_P6"
gene_integrated$phate_cluster[which(gene_integrated$new_ident %in% "CTCL_Patients_2")] <- "CTCL_P11"
gene_integrated$phate_cluster[which(gene_integrated$new_ident %in% "CTCL_Patients_3")] <- "CTCL_P3"
gene_integrated$phate_cluster[which(gene_integrated$new_ident %in% "CTCL_Patients_4")] <- "CTCL_P1"
gene_integrated$phate_cluster[which(gene_integrated$new_ident %in% "CTCL_Patients_5")] <- "CTCL_P4"
gene_integrated$phate_cluster[which(gene_integrated$new_ident %in% "CTCL_Patients_6")] <- "CTCL_P8"
gene_integrated$phate_cluster[which(gene_integrated$new_ident %in% "CTCL_Patients_7")] <- "CTCL_P9"
gene_integrated$phate_cluster[which(gene_integrated$new_ident %in% "CTCL_Patients_8")] <- "CTCL_P2"
gene_integrated$phate_cluster[which(gene_integrated$new_ident %in% "CTCL_Patients_9")] <- "CTCL_P7"
gene_integrated$phate_cluster[which(gene_integrated$new_ident %in% "CTCL_Patients_10")] <- "CTCL_P10"
gene_integrated$phate_cluster[which(gene_integrated$new_ident %in% "CTCL_Patients_11")] <- "CTCL_P5"

Three_group_cluster <- gene_integrated
Three_group_cluster$three_clusters <- "Group1"
Three_group_cluster$three_clusters[which(Three_group_cluster$phate_cluster %in% "CTCL_P5")] <- "Group3"
Three_group_cluster$three_clusters[which(Three_group_cluster$phate_cluster %in% "CTCL_P10")] <- "Group3"
Three_group_cluster$three_clusters[which(Three_group_cluster$phate_cluster %in% "CTCL_P7")] <- "Group3"
Three_group_cluster$three_clusters[which(Three_group_cluster$phate_cluster %in% "CTCL_P9")] <- "Group3"
Three_group_cluster$three_clusters[which(Three_group_cluster$phate_cluster %in% "CTCL_P2")] <- "Group3"
Three_group_cluster$three_clusters[which(Three_group_cluster$phate_cluster %in% "CTCL_P8")] <- "Group3"
Three_group_cluster$three_clusters[which(Three_group_cluster$phate_cluster %in% "CTCL_P6")] <- "Group3"
Three_group_cluster$three_clusters[which(Three_group_cluster$phate_cluster %in% "CTCL_P11")] <- "Group3"
Three_group_cluster$three_clusters[which(Three_group_cluster$phate_cluster %in% "CTCL_P3")] <- "Group3"
Three_group_cluster$three_clusters[which(Three_group_cluster$phate_cluster %in% "CTCL_P1")] <- "Group3"
Three_group_cluster$three_clusters[which(Three_group_cluster$phate_cluster %in% "CTCL_P4")] <- "Group3"
Three_group_cluster$three_clusters[which(Three_group_cluster$phate_cluster %in% "CTCL_Normal")] <- "Group2"

Idents(Three_group_cluster) <- "three_clusters"
saveRDS(Three_group_cluster, file = file.path("~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/UMAP/", "gene_integrated_Three_group_cluster.rds"))

#Figure3A
Three_group_cluster <- readRDS(file ="~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/UMAP/gene_integrated_Three_group_cluster.rds")
DimPlot(Three_group_cluster, reduction = "umap", label = FALSE,group.by = "three_clusters", split.by = "category",
        cols = c("Group1" = "#1E90FF", "Group2" = "#00BFFF", "Group3" = "#FF0000"))##, group.by = "three_clusters", split.by = "category" )


Idents(gene_integrated) <- "phate_cluster"
a <- subset(gene_integrated, ident = "CTCL_Normal")
gene_integrated$phate_cluster[which(a$clonotype_id == "clonotype1")] <- "CTCL_Nomal_clonotype1"
barcode_clonotype1_in_normal <- colnames(gene_integrated)[which(a$clonotype_id == "clonotype1")]
gene_integrated$phate_cluster[which(colnames(gene_integrated) %in% barcode_clonotype1_in_normal)] <- "CTCL_Nomal_clonotype1"

gene_integrated_with_transit <- gene_integrated 
saveRDS(gene_integrated_with_transit, file = file.path("~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/UMAP/", "gene_integrated_with_transit.rds"))

gene_integrated_with_transit <- readRDS("~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/UMAP/gene_integrated_with_transit.rds")
DimPlot(gene_integrated_with_transit, reduction = "umap", split.by = "category", label = TRUE, group.by = "phate_cluster")
Idents(gene_integrated_with_transit) <- "category"
CTCL_patient <- subset(gene_integrated_with_transit, ident = "CTCL_Patients")
##separate the integrated CTCL_patient into three cluster
CTCL_patient$phate_cluster_CTCL_patient <- "Normal_CD4T normal_clonotype1-"## "Normal_CD4T normal_clonotype1+"  "CTCL CTCL"
CTCL_cluster <- c("CTCL_P8", "CTCL_P7", "CTCL_P6", "CTCL_P4", "CTCL_P3", "CTCL_P2", "CTCL_P1","CTCL_P11", "CTCL_P10", "CTCL_P9", "CTCL_P5")
CTCL_patient$phate_cluster_CTCL_patient[which(CTCL_patient$phate_cluster %in% CTCL_cluster)] <- "CTCL CTCL"
CTCL_patient$phate_cluster_CTCL_patient[which(CTCL_patient$phate_cluster == "CTCL_Nomal_clonotype1")] <- "Normal_CD4T normal_clonotype1+"
head(CTCL_patient@meta.data)
Idents(CTCL_patient) <- "phate_cluster_CTCL_patient"

DimPlot(CTCL_patient, reduction = "umap",group.by = "phate_cluster_CTCL_patient")#raster=FALSE

saveRDS(CTCL_patient, file = file.path("~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/UMAP/", "CTCL_patient_integrated.rds"))
#Figure 1C
CTCL_patient <- readRDS("~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/UMAP/CTCL_patient_integrated.rds")
Idents(CTCL_patient) <- "seurat_clusters"
cols_mix <- c("0"= "#00BFFF", "1"= "#D8BFD8", "2"= "#FF6347", "3"= "#FFC0CB" ,"4"= "#FF4500", "5" = "#FF8C00", "6"= "#DAA520", "7"="#6A5ACD", "8"="#8B008B", "9"="#FF00FF", "10"= "#FFD700", "11"= "#DA70D6", "12"= "#00BFFF", "13"= "#00BFFF")
DimPlot(CTCL_patient, reduction = "umap", group.by = "seurat_clusters", cols = cols_mix, label = FALSE) + NoLegend()

#Figure 1D
CTCL_patient$P <- "others"
CTCL_patient$P[which(CTCL_patient$orig.ident == "P_5P_tcr")] <- "P_others"
CTCL_patient$P[which(CTCL_patient$orig.ident == "P_5P_tcr" & CTCL_patient$clonotype_id %in% c("clonotype1", "clonotype2", "clonotype3", "clonotype4"))] <- "P7_CTCL"

clonotype1 <- subset (CTCL_patient, subset = P == "P_CTCL")
clonotype2 <- subset (CTCL_patient, subset = P == "P_others")
clonotype1 <- Cells(clonotype1)
clonotype2 <- Cells(clonotype2)
DimPlot (CTCL_patient, cells.highlight = list(clonotype1, clonotype2), cols.highlight = c("#00BFFF","#FF0000"), reduction = "umap")+NoLegend()  

#Figure 1E
#Th2-like: GATA3>1, CCR4>1; Treg-like: FOXP3>1, CTLA4>1; Tfh-like: PDCD1>1, CXCR5>1; Th17-like: RORC>1, CCR6>1;Th1-like:TBX21>1 & IFNG >1
DefaultAssay(CTCL_patient) <- "RNA"
Th2-like <- subset(CTCL_patient, subset = GATA3 >1 & CCR4 >1 )
Th2-like <- Cells(Th2-like)
Th2-like <- DimPlot(CTCL_patient, cells.highlight = Th2-like, reduction = "umap",label = FALSE,pt.size = 0.01)+NoLegend()
Th2-like
Treg-like <- subset(CTCL_patient, subset = FOXP3 >1 & CTLA4 >1 )
Treg-like <- Cells(Treg-like)
Treg-like <- DimPlot(CTCL_patient, cells.highlight = Treg-like, reduction = "umap",label = FALSE,pt.size = 0.01)+NoLegend()
Treg-like
Tfh-like <- subset(CTCL_patient, subset = PDCD1 >1 & CXCR5 >1 )
Tfh-like <- Cells(Tfh-like)
Tfh-like <- DimPlot(CTCL_patient, cells.highlight = Tfh-like, reduction = "umap",label = FALSE,pt.size = 0.01)+NoLegend()
Tfh-like
Th17-like <- subset(CTCL_patient, subset = RORC>1 >1 & CCR6 >1 )
Th17-like <- Cells(Th17-like)
Th17-like <- DimPlot(CTCL_patient, cells.highlight = Th17-like, reduction = "umap",label = FALSE,pt.size = 0.01)+NoLegend()
Th17-like
Th1-like <- subset(CTCL_patient, subset = TBX21>1 & IFNG >1 )
Th1-like <- Cells(Th1-like)
Th1-like <- DimPlot(CTCL_patient, cells.highlight = Th1-like, reduction = "umap",label = FALSE,pt.size = 0.01)+NoLegend()
Th1-like

#Figure 4A
Three_group_cluster <- readRDS(file ="~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/UMAP/gene_integrated_Three_group_cluster.rds")
DimPlot(Three_group_cluster, reduction = "umap", label = FALSE,group.by = "three_clusters", split.by = "category",
        cols = c("Group1" = "#1E90FF", "Group2" = "#00BFFF", "Group3" = "#FF0000"))##, group.by = "three_clusters", split.by = "category" )
Idents(Three_group_cluster) <- "three_clusters"
DefaultAssay(Three_group_cluster) <- "RNA"
CTCL_NP.NH  <- FindMarkers(Three_group_cluster, ident.1 = "Group3", ident.2 = c("Group1", "Group2"))#Compare CTCL clusters to both normal CD4T cells in the patients(NP) and CD4 T cells in healthy controls(NH): CTCL_NP.NH 
write.csv(CTCL_NP.NH, file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/Result/DEG/CTCL_NP.NH.csv")
##Heatmap
data.S <- Three_group_cluster
DefaultAssay(data.S)<- "RNA"
data.S <- FindVariableFeatures(data.S, nfeatures =19371)  #in order to plot all of genes in the Doheatmap
data.S <- ScaleData(data.S)
features <- read.csv(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/Result/DEG/CTCL_NP.NH_filtered.csv", row.names = 1)#heatmap: min. pct > 0.25, min. diff. pct > 0.2, adjusted p value d 0.05, logfc threshold = 0.25
DoHeatmap(subset(data.S, downsample = 10000),slot = "scale.data",features = as.character(features$X), size = 3, group.by =  "three_clusters", label = FALSE)+scale_fill_gradientn(colors = c("blue", "white", "red"))+NoLegend()+ theme(text = element_text(size = 0))
subset(data.S)
DoHeatmap(data.S,slot = "scale.data",features = as.character(features$X), size = 3, group.by =  "three_clusters", label = FALSE)+scale_fill_gradientn(colors = c("blue", "white", "red"))+NoLegend()+ theme(text = element_text(size = 0))
tiff("~/Desktop/heatmap.tiff", units="in", width=4, height=7, res=300)
print(DoHeatmap(subset(data.S, downsample = 10000),slot = "scale.data",features = as.character(features$X), size = 3, group.by =  "three_clusters", label = FALSE)+scale_fill_gradientn(colors = c("blue", "white", "red"))+NoLegend()+ theme(text = element_text(size = 6)))
dev.off()
#Figure 4C
VlnPlot(data.S, features = "CCR7", pt.size = 0, cols =c("blue", "red"))+NoLegend()+scale_y_continuous(limits = c(0,3))


#Supplementary Figure
CTCLs_11_DEG <- FindAllMarkers(CTCLs_clusters)
write.csv(CTCLs_11_DEG, file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/Result/DEG/CTCLs_11_DEG.csv")
## plot: heatmap
data.S <- CTCLs_clusters
DefaultAssay(data.S)<- "RNA"
data.S <- FindVariableFeatures(data.S, nfeatures =19371)  #in order to plot all of genes in the Doheatmap
data.S <- ScaleData(data.S)
features <- read.csv(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/Result/DEG/CTCLs_11_DEG.csv", row.names = 1)
pos_features <- filter(features, features$avg_logFC > 0)
pos_top10_features <- pos_features %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_logFC)
DoHeatmap(subset(data.S, downsample = 1000),slot = "scale.data",features = as.character(pos_top10_features$gene), size = 3, group.by =  "seurat_clusters")+scale_fill_gradientn(colors = c("blue", "white", "red"))














##Phate only for CTCL_patients
#1. take clonotype1+ cells in the normal CD4 T cells as a transitional cluster
CTCL_patient <- readRDS("~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/UMAP/CTCL_patient_integrated.rds")
##single patient annotation
Idents(CTCL_patient) <- "seurat_clusters"
cols_mix <- c("0"= "#00BFFF", "1"= "#D8BFD8", "2"= "#FF6347", "3"= "#FFC0CB" ,"4"= "#FF4500", "5" = "#FF8C00", "6"= "#DAA520", "7"="#6A5ACD", "8"="#8B008B", "9"="#FF00FF", "10"= "#FFD700", "11"= "#DA70D6", "12"= "#00BFFF", "13"= "#00BFFF")
DimPlot(CTCL_patient, reduction = "umap", group.by = "seurat_clusters", cols = cols_mix, label = FALSE) + NoLegend()
data.sample <- subset(CTCL_patient, subset = orig.ident == "SEOG_5P_tcr")#based on patient_ID
data.sample <- Cells(data.sample) 
DimPlot(CTCL_patient, cells.highlight = data.sample, reduction = "umap", label = FALSE,pt.size = 0.02 )+NoLegend()

data.sample <- subset(CTCL_patient, subset = cdr3s_aa == "TRA:CAVSEPGGYQKVTF;TRB:CASRAQPGLAGNEQFF")#based on top clonotype cdr3_aa sequence
data.sample <- Cells(data.sample) 
DimPlot(CTCL_patient, cells.highlight = data.sample, reduction = "umap", label = FALSE )+NoLegend() #,pt.size = 0.02

###Figure1
sample <- subset(CTCL_patient, subset = orig.ident == "NOV1819_5P_tcr")
sample <- Cells(sample)
sample <- DimPlot(CTCL_patient, cells.highlight = sample, reduction = "umap", label = FALSE ) ##+NoLegend()
sample
sample <- DimPlot(CTCL_patient, cells.highlight = sample, reduction = "umap", label = FALSE ,pt.size = 0.02, cols.highlight = "#FF6347",cols="#A9A9A9" )  ##+ scale_color_manual(labels = c("other", "sample"), values = c("#A9A9A9", "#FF6347")) 
sample

treg <- subset(CTCL_patient, subset = RORC >1 & CCR6 >1 )
treg <- Cells(treg)
treg_plot <- DimPlot(CTCL_patient, cells.highlight = treg, reduction = "umap",label = FALSE,pt.size = 0.01)+NoLegend()
treg_plot

##DEG
Idents(CTCL_patient) <- "phate_cluster_CTCL_patient"
DefaultAssay(CTCL_patient) <- "RNA"

##Find monototic trend among transitional cluster
CTCL_patient_C1 <- FindMarkers(CTCL_patient, ident.1 = "CTCL CTCL", ident.2 = "Normal_CD4T normal_clonotype1+")
CTCL_patient_C2 <- FindMarkers(CTCL_patient, ident.1 = "CTCL CTCL", ident.2 = "Normal_CD4T normal_clonotype1-")
CTCL_patient_C3 <- FindMarkers(CTCL_patient, ident.1 = "Normal_CD4T normal_clonotype1+", ident.2 = "Normal_CD4T normal_clonotype1-")
write.csv(CTCL_patient_C1, file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/Result/DEG/CTCL_patient_C1.csv")
write.csv(CTCL_patient_C2, file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/Result/DEG/CTCL_patient_C2.csv")
write.csv(CTCL_patient_C3, file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/Result/DEG/CTCL_patient_C3.csv")
##plot heatmap

data.S <- CTCL_patient
DefaultAssay(data.S)<- "RNA"
data.S <- FindVariableFeatures(data.S, nfeatures =19371)  #in order to plot all of genes in the Doheatmap
data.S <- ScaleData(data.S)
Idents(data.S) <- "phate_cluster_CTCL_patient"
features <- read.csv(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/Result/DEG/monototical_features_CTCL_patients.csv", row.names = 1)
DoHeatmap(subset(data.S, downsample = 10000),slot = "scale.data",features = as.character(features$X), size = 3, group.by =  "phate_cluster_CTCL_patient")
DoHeatmap(subset(data.S, downsample = 10000),slot = "scale.data",features = as.character(features$x), size = 3, group.by =  "phate_cluster_CTCL_patient")+scale_fill_gradientn(colors = c("blue", "white", "red"))

cluster.averages <- AverageExpression(data.S, return.seurat = TRUE) ##plot average expression of genes
DoHeatmap(cluster.averages,features = as.character(features$x),size = 3, draw.lines = FALSE)+scale_fill_gradientn(colors = c("blue", "white", "red"))   #modify the size of the text: +theme(text = element_text(size = 20))



#2. take cells identified as transitional cells in each single patient analysis as a transitional cluster
CTCL_patient_2 <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/UMAP/CTCL_patient_2_integrated.rds")
col_three_cluster <- c("CTCL CTCL" = "#FF0000", "Normal_CD4T normal_clonotype1-"= "#0000FF", "Normal_CD4T normal_clonotype1+" = "#8B008B")
DimPlot(CTCL_patient_2, reduction = "umap",group.by = "phate_cluster_CTCL_patient", label = FALSE, cols = c("CTCL CTCL" = "#FF0000", "Normal_CD4T normal_clonotype1-"= "#0000FF", "Normal_CD4T normal_clonotype1+" = "#8B008B"))+NoLegend()
Idents(CTCL_patient_2) <- "phate_cluster_CTCL_patient"
DefaultAssay(CTCL_patient_2) <- "RNA"
CTCL_pateint_2_C1 <- FindMarkers(CTCL_patient_2, ident.1 = "CTCL CTCL", ident.2 =  "Normal_CD4T normal_clonotype1-")
CTCL_pateint_2_C2 <- FindMarkers(CTCL_patient_2, ident.1 = "CTCL CTCL", ident.2 =  "Normal_CD4T normal_clonotype1+")
CTCL_pateint_2_C3 <- FindMarkers(CTCL_patient_2, ident.1 = "Normal_CD4T normal_clonotype1+", ident.2 = "Normal_CD4T normal_clonotype1-")
write.csv(CTCL_pateint_2_C1, file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/Result/DEG/CTCL_pateint_2_C1.csv")
write.csv(CTCL_pateint_2_C2, file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/Result/DEG/CTCL_pateint_2_C2.csv")
write.csv(CTCL_pateint_2_C3, file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/Result/DEG/CTCL_pateint_2_C3.csv")
##plot heatmap
data.S <- CTCL_patient_2
DefaultAssay(data.S)<- "RNA"
data.S <- FindVariableFeatures(data.S, nfeatures =19371)  #in order to plot all of genes in the Doheatmap
data.S <- ScaleData(data.S)
Idents(data.S) <- "phate_cluster_CTCL_patient"

VlnPlot(data.S, features = c("CD82", "CCR4", "KIR3DL2"), pt.size = 0)+NoLegend() ##Biomarker CD82, CCR4. KIR3DL2
DotPlot(data.S, features = c("CD82", "CCR4", "KIR3DL2"), cols = c("blue", "red"))+RotatedAxis()+NoLegend()



features <- read.csv(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/Result/DEG/monototical_features_CTCL_patients_2.csv", row.names = 1)
DoHeatmap(subset(data.S, downsample = 10000),slot = "scale.data",features = as.character(features$x), size = 3, group.by =  "phate_cluster_CTCL_patient")
DoHeatmap(subset(data.S, downsample = 10000),slot = "scale.data",features = c("CCR4", "CD82",as.character(features$x)), size = 3, group.by =  "phate_cluster_CTCL_patient")+scale_fill_gradientn(colors = c("blue", "white", "red"))

cluster.averages <- AverageExpression(data.S, return.seurat = TRUE) ##plot average expression of genes
DoHeatmap(cluster.averages,features = c("CCR4", "CD82", as.character(features$x) ),size = 3, draw.lines = FALSE)+scale_fill_gradientn(colors = c("blue", "white", "red"))   #modify the size of the text: +theme(text = element_text(size = 20))
##*****use the DEG of CTCL_NP.NH_filtered to plot******######
features <- read.csv(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/Result/DEG/CTCL_NP.NH_filtered.csv", row.names = 1)
DoHeatmap(subset(data.S, downsample = 10000),slot = "scale.data",features = as.character(features$X), size = 3, group.by =  "phate_cluster_CTCL_patient")+scale_fill_gradientn(colors = c("blue", "white", "red"))

cluster.averages <- AverageExpression(data.S, return.seurat = TRUE)
DoHeatmap(cluster.averages,features = as.character(features$X), size = 3, draw.lines = FALSE)+scale_fill_gradientn(colors = c("blue", "white", "red"))+NoLegend()
tiff("~/Desktop/heatmap_1.tiff", units="in", width=4, height=7, res=300)
print(DoHeatmap(cluster.averages,features = as.character(features$X), size = 3, draw.lines = FALSE)+scale_fill_gradientn(colors = c("blue", "white", "red"))+NoLegend()+theme(text = element_text(size = 6)))
dev.off()

##DEG for 11 different patients CTCL cells
CTCLs_clusters <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/UMAP/CTCLs_clusters.rds")
DimPlot(CTCLs_clusters, reduction = "umap", split.by = "category", label = TRUE )



##DEG1: SS vs MF
DefaultAssay(CTCLs_clusters)<- "RNA"
CTCL_cluster_DEG1 <- FindMarkers(CTCLs_clusters, ident.1 = c(0,2,3,4,5,8,10), ident.2 = c(1,6,7,9))
write.csv(CTCL_cluster_DEG1, file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/Result/DEG/CTCL_cluster_DEG1.csv")
##DEG2: ISCL stage III or more vs II or less 
CTCL_cluster_DEG2 <- FindMarkers(CTCLs_clusters, ident.1 = c(0,2,3,4,5,8,10), ident.2 = c(1,6,7,9))
##DEG3: B stage, B2 vs B1
CTCL_cluster_DEG3 <- FindMarkers(CTCLs_clusters, ident.1 = c(0,2,3,4,8,10), ident.2 = c(1,5,6,7,9))
write.csv(CTCL_cluster_DEG3, file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/Result/DEG/CTCL_cluster_DEG3.csv")
CTCLs_clusters$DEG1 <- "SS"
CTCLs_clusters$DEG1[which(CTCLs_clusters$seurat_clusters %in% c(1,6,7,9))] <- "MF"
##ploting DEG1&2 and DEG3
#####create a new column in the metadata
CTCLs_clusters$DEG1 <- "SS"
CTCLs_clusters$DEG1[which(CTCLs_clusters$seurat_clusters %in% c(1,6,7,9))] <- "MF"
CTCLs_clusters$DEG3 <- "B2"
CTCLs_clusters$DEG3[which(CTCLs_clusters$seurat_clusters %in% c(1,5, 6,7,9))] <- "B1"
Idents(CTCLs_clusters) <- "DEG1"
Idents(CTCLs_clusters) <- "DEG3"
data.S <- CTCLs_clusters
DefaultAssay(data.S)<- "RNA"
data.S <- FindVariableFeatures(data.S, nfeatures =19371)  #in order to plot all of genes in the Doheatmap
data.S <- ScaleData(data.S)

features <- read.csv(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/Result/DEG/CTCL_cluster_DEG1_filtered.csv", row.names = 1)
Idents(data.S) <- "DEG1"
DoHeatmap(subset(data.S, downsample = 10000),slot = "scale.data",features = as.character(features$X), size = 3, group.by =  "DEG1")+scale_fill_gradientn(colors = c("blue", "white", "red"))
cluster.averages <- AverageExpression(data.S, return.seurat = TRUE) ##plot average expression of genes
DoHeatmap(cluster.averages, features = as.character(features$X),draw.lines = FALSE, size = 3)+scale_fill_gradientn(colors = c("blue", "white", "red"))  #modify the size of the text: +theme(text = element_text(size = 20))

StackedVlnPlot(data.S, features = as.character(features$X))

tiff("~/Desktop/heatmap1.tiff", units="in", width=3, height=5, res=300)
print(DoHeatmap(subset(data.S, downsample = 10000),slot = "scale.data",features = as.character(features$X), size = 3, group.by =  "DEG1")+scale_fill_gradientn(colors = c("blue", "white", "red"))+NoLegend()+ theme(text = element_text(size = 6)))
dev.off()
tiff("~/Desktop/heatmap2.tiff", units="in", width=3, height=5, res=300)
print(DoHeatmap(cluster.averages, features = as.character(features$X),draw.lines = FALSE, size = 3)+scale_fill_gradientn(colors = c("blue", "white", "red")))  #modify the size of the text: +theme(text = element_text(size = 20))
dev.off()
tiff("~/Desktop/vlnplot.tiff", units="in", width=5, height=10, res=300)
print(StackedVlnPlot(data.S, features = as.character(features$X)))
dev.off()






CTCLs_11_DEG <- FindAllMarkers(CTCLs_clusters)
write.csv(CTCLs_11_DEG, file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/Result/DEG/CTCLs_11_DEG.csv")
## plot: heatmap
data.S <- CTCLs_clusters
DefaultAssay(data.S)<- "RNA"
data.S <- FindVariableFeatures(data.S, nfeatures =19371)  #in order to plot all of genes in the Doheatmap
data.S <- ScaleData(data.S)
features <- read.csv(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/Result/DEG/CTCLs_11_DEG.csv", row.names = 1)
DoHeatmap(subset(data.S, downsample = 10000),slot = "scale.data",features = as.character(features$X), size = 3, group.by =  "seurat_clusters")
DoHeatmap(subset(data.S, downsample = 1000),slot = "scale.data",features = as.character(rownames(features)), size = 3, group.by =  "seurat_clusters")+scale_fill_gradientn(colors = c("blue", "white", "red"))
pos_features <- filter(features, features$avg_logFC > 0)
pos_top10_features <- pos_features %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_logFC)
DoHeatmap(subset(data.S, downsample = 1000),slot = "scale.data",features = as.character(pos_top10_features$gene), size = 3, group.by =  "seurat_clusters")+scale_fill_gradientn(colors = c("blue", "white", "red"))


##DEG for differential categoriezed by severity of 11 group of CTCL cells: S vs M+L 
S_M.L_markers <- FindMarkers(CTCLs_clusters, ident.1 = c(0,3,4,5,8,10), ident.2 = c(1,2,6,7,9))
write.csv(CTCLs_S_M.L_markers, file ="CTCLs_~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/Result/DEG/CTCLs_S_M.L_markers.csv")
data.S$severity <- "S"
data.S$severity[which(data.S$seurat_clusters %in% c(1,2,9))] <- "M"
data.S$severity[which(data.S$seurat_clusters %in% c(6,7))] <- "L"
Idents(data.S) <- "severity"
###plot heatmap
features <- read.csv(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/Result/DEG/CTCLs_S_M.L_markers_filtered.csv", row.names = 1)
DoHeatmap(subset(data.S, downsample = 5000),slot = "scale.data",features = as.character(features$X), size = 3, group.by = "severity")+scale_fill_gradientn(colors = c("blue", "white", "red"))

cluster.averages <- AverageExpression(data.S, return.seurat = TRUE) ##plot average expression of genes
DoHeatmap(cluster.averages,features = as.character(features$X),size = 3, draw.lines = FALSE)+scale_fill_gradientn(colors = c("blue", "white", "red"))   #modify the size of the text: +theme(text = element_text(size = 20))

DotPlot(data.S, features = as.character(features$X), cols = c("#FFD700", "#008000"))+ RotatedAxis() 


CTCL_patient <- readRDS("~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Integration/UMAP/CTCL_patient_integrated.rds")


CTCL_patient$P7 <- "others"
CTCL_patient$P7[which(CTCL_patient$orig.ident == "FEB112020_2_5P_tcr")] <- "P7_others"
CTCL_patient$P7[which(CTCL_patient$orig.ident == "FEB112020_2_5P_tcr" & CTCL_patient$clonotype_id %in% c("clonotype1", "clonotype2", "clonotype3", "clonotype4"))] <- "P7_CTCL"

clonotype1 <- subset (CTCL_patient, subset = P7 == "P7_CTCL")
clonotype2 <- subset (CTCL_patient, subset = P7 == "P7_others")
clonotype1 <- Cells(clonotype1)
clonotype2 <- Cells(clonotype2)
DimPlot (CTCL_patient, cells.highlight = list(clonotype1, clonotype2), cols.highlight = c("#00BFFF","#FF0000"), reduction = "umap")+NoLegend()  ##+scale_color_manual(labels = c("other", "clonotype1"), values = c("grey", "red"))
#SEOG
CTCL_patient$P1 <- "others"
CTCL_patient$P1[which(CTCL_patient$orig.ident == "SEOG_5P_tcr")] <- "P1_others"
CTCL_patient$P1[which(CTCL_patient$orig.ident == "SEOG_5P_tcr" & CTCL_patient$clonotype_id %in% c("clonotype1", "clonotype2", "clonotype3"))] <- "P1_CTCL"

clonotype1 <- subset (CTCL_patient, subset = P1 == "P1_CTCL")
clonotype2 <- subset (CTCL_patient, subset = P1 == "P1_others")
clonotype1 <- Cells(clonotype1)
clonotype2 <- Cells(clonotype2)
DimPlot (CTCL_patient, cells.highlight = list(clonotype1, clonotype2), cols.highlight = c("#00BFFF","#FF0000"), reduction = "umap")+NoLegend()  ##+scale_color_manual(labels = c("other", "clonotype1"), values = c("grey", "red"))
#JA1819
CTCL_patient$P2 <- "others"
CTCL_patient$P2[which(CTCL_patient$orig.ident == "JA1819_5P_tcr")] <- "P2_others"
CTCL_patient$P2[which(CTCL_patient$orig.ident == "JA1819_5P_tcr" & CTCL_patient$clonotype_id %in% c("clonotype1", "clonotype2", "clonotype3"))] <- "P2_CTCL"

clonotype1 <- subset (CTCL_patient, subset = P2 == "P2_CTCL")
clonotype2 <- subset (CTCL_patient, subset = P2 == "P2_others")
clonotype1 <- Cells(clonotype1)
clonotype2 <- Cells(clonotype2)
DimPlot (CTCL_patient, cells.highlight = list(clonotype1, clonotype2), cols.highlight = c("#00BFFF","#FF0000"), reduction = "umap")+NoLegend()  ##+scale_color_manual(labels = c("other", "clonotype1"), values = c("grey", "red"))
#NOV1819
CTCL_patient$P3 <- "others"
CTCL_patient$P3[which(CTCL_patient$orig.ident == "NOV1819_5P_tcr")] <- "P3_others"
CTCL_patient$P3[which(CTCL_patient$orig.ident == "NOV1819_5P_tcr" & CTCL_patient$clonotype_id %in% c("clonotype1", "clonotype2"))] <- "P3_CTCL"

clonotype1 <- subset (CTCL_patient, subset = P3 == "P3_CTCL")
clonotype2 <- subset (CTCL_patient, subset = P3 == "P3_others")
clonotype1 <- Cells(clonotype1)
clonotype2 <- Cells(clonotype2)
DimPlot (CTCL_patient, cells.highlight = list(clonotype1, clonotype2), cols.highlight = c("#00BFFF","#FF0000"), reduction = "umap")+NoLegend()  ##+scale_color_manual(labels = c("other", "clonotype1"), values = c("grey", "red"))

#A23
CTCL_patient$P4 <- "others"
CTCL_patient$P4[which(CTCL_patient$orig.ident == "A23_5P_tcr")] <- "P4_others"
CTCL_patient$P4[which(CTCL_patient$orig.ident == "A23_5P_tcr" & CTCL_patient$clonotype_id %in% c("clonotype1", "clonotype2", "clonotype3"))] <- "P4_CTCL"

clonotype1 <- subset (CTCL_patient, subset = P4 == "P4_CTCL")
clonotype2 <- subset (CTCL_patient, subset = P4 == "P4_others")
clonotype1 <- Cells(clonotype1)
clonotype2 <- Cells(clonotype2)
DimPlot (CTCL_patient, cells.highlight = list(clonotype1, clonotype2), cols.highlight = c("#00BFFF","#FF0000"), reduction = "umap")+NoLegend()  ##+scale_color_manual(labels = c("other", "clonotype1"), values = c("grey", "red"))
#A14
CTCL_patient$P5 <- "others"
CTCL_patient$P5[which(CTCL_patient$orig.ident == "A14_5P_tcr")] <- "P5_others"
CTCL_patient$P5[which(CTCL_patient$orig.ident == "A14_5P_tcr" & CTCL_patient$clonotype_id %in% c("clonotype1", "clonotype2", "clonotype3"))] <- "P5_CTCL"

clonotype1 <- subset (CTCL_patient, subset = P5 == "P5_CTCL")
clonotype2 <- subset (CTCL_patient, subset = P5 == "P5_others")
clonotype1 <- Cells(clonotype1)
clonotype2 <- Cells(clonotype2)
DimPlot (CTCL_patient, cells.highlight = list(clonotype1, clonotype2), cols.highlight = c("#00BFFF","#FF0000"), reduction = "umap")+NoLegend()  ##+scale_color_manual(labels = c("other", "clonotype1"), values = c("grey", "red"))

#MA1319
CTCL_patient$P6 <- "others"
CTCL_patient$P6[which(CTCL_patient$orig.ident == "MA1319_5P_tcr")] <- "P6_others"
CTCL_patient$P6[which(CTCL_patient$orig.ident == "MA1319_5P_tcr" & CTCL_patient$clonotype_id %in% c("clonotype1", "clonotype2", "clonotype3"))] <- "P6_CTCL"

clonotype1 <- subset (CTCL_patient, subset = P6 == "P6_CTCL")
clonotype2 <- subset (CTCL_patient, subset = P6 == "P6_others")
clonotype1 <- Cells(clonotype1)
clonotype2 <- Cells(clonotype2)
DimPlot (CTCL_patient, cells.highlight = list(clonotype1, clonotype2), cols.highlight = c("#00BFFF","#FF0000"), reduction = "umap")+NoLegend()  ##+scale_color_manual(labels = c("other", "clonotype1"), values = c("grey", "red"))

#FEB192020
CTCL_patient$P8 <- "others"
CTCL_patient$P8[which(CTCL_patient$orig.ident == "FEB192020_5P_tcr")] <- "P8_others"
CTCL_patient$P8[which(CTCL_patient$orig.ident == "FEB192020_5P_tcr" & CTCL_patient$clonotype_id %in% c("clonotype1", "clonotype2", "clonotype3"))] <- "P8_CTCL"

clonotype1 <- subset (CTCL_patient, subset = P8 == "P8_CTCL")
clonotype2 <- subset (CTCL_patient, subset = P8 == "P8_others")
clonotype1 <- Cells(clonotype1)
clonotype2 <- Cells(clonotype2)
DimPlot (CTCL_patient, cells.highlight = list(clonotype1, clonotype2), cols.highlight = c("#00BFFF","#FF0000"), reduction = "umap")+NoLegend()  ##+scale_color_manual(labels = c("other", "clonotype1"), values = c("grey", "red"))

#JUL18
CTCL_patient$P9 <- "others"
CTCL_patient$P9[which(CTCL_patient$orig.ident == "JUL18_5P_tcr")] <- "P9_others"
CTCL_patient$P9[which(CTCL_patient$orig.ident == "JUL18_5P_tcr" & CTCL_patient$clonotype_id %in% c("clonotype1", "clonotype2", "clonotype3"))] <- "P9_CTCL"

clonotype1 <- subset (CTCL_patient, subset = P9 == "P9_CTCL")
clonotype2 <- subset (CTCL_patient, subset = P9 == "P9_others")
clonotype1 <- Cells(clonotype1)
clonotype2 <- Cells(clonotype2)
DimPlot (CTCL_patient, cells.highlight = list(clonotype1, clonotype2), cols.highlight = c("#00BFFF","#FF0000"), reduction = "umap")+NoLegend()  ##+scale_color_manual(labels = c("other", "clonotype1"), values = c("grey", "red"))
#AUG1914
CTCL_patient$P11 <- "others"
CTCL_patient$P11[which(CTCL_patient$orig.ident == "AUG1914_5P_tcr")] <- "P11_others"
CTCL_patient$P11[which(CTCL_patient$orig.ident == "AUG1914_5P_tcr" & CTCL_patient$clonotype_id %in% c("clonotype1", "clonotype2", "clonotype3", "clonotype4", "clonotype5", "clonotype6", "clonotype7"))] <- "P11_CTCL"

clonotype1 <- subset (CTCL_patient, subset = P11 == "P11_CTCL")
clonotype2 <- subset (CTCL_patient, subset = P11 == "P11_others")
clonotype1 <- Cells(clonotype1)
clonotype2 <- Cells(clonotype2)
DimPlot (CTCL_patient, cells.highlight = list(clonotype1, clonotype2), cols.highlight = c("#00BFFF","#FF0000"), reduction = "umap")+NoLegend()  ##+scale_color_manual(labels = c("other", "clonotype1"), values = c("grey", "red"))


##Figure 4C
##InferCNV
P_phate <- readRDS(file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/Phate_subpopulations/rsd_file/P_phate.rds" )
DimPlot(P_phate, reduction = "tsne", pt.size = 0.5)
write.table(as.matrix(GetAssayData(object = P_phate, slot = "counts")), 
            '~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/20210621/seurat4-202101011/inferCNV_20211020/P/P_phate_counts_matrix.txt', 
            sep = '\t', row.names = T, col.names = T, quote = F)
P_phate$inferCNV <- "CTCL"
P_phate$inferCNV[which(P_phate$subseted_cluster %in% "Normal_CD4T normal_clonotype1-")] <- "Normal"
P_phate$inferCNV[which(P_phate$subseted_cluster %in% "Normal_CD4T normal_clonotype1+")] <- "Transitory"
write.table(P_phate$inferCNV, file = "~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/20210621/seurat4-202101011/inferCNV_20211020/P/P_phate_metadata.txt", sep = '\t', row.names = T, col.names = F, quote = F)


library(infercnv)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/20210621/seurat4-202101011/inferCNV_20211020/P/P_phate_counts_matrix.txt",
                                    annotations_file="~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/20210621/seurat4-202101011/inferCNV_20211020/P/P_phate_metadata.txt",
                                    delim="\t",
                                    gene_order_file= "./refdata-cellranger-GRCh38-3.0.0_gene_order_v4.txt",
                                    ref_group_names=c("Normal" )) 
infercnv_obj_3 = infercnv::run(infercnv_obj,
                               cutoff=0.1, 
                               out_dir="~/project/2021_CTCL_data/CTCL_subpopolation_investigation/20210329/20210621/seurat4-202101011/inferCNV_20211020/P", 
                               cluster_by_groups=T,   
                               denoise=T,
                               HMM=T,
                               plot_steps=FALSE,
                               no_prelim_plot=TRUE,
                               png_res=60,
                               HMM_type='i3')

