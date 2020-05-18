library(Seurat)
library(hdf5r)
library(sctransform)
library(ggplot2)
setwd("/home/rstudio/scRNAseq_COVID/")

####set seed
set.seed(2020)

## download count matrices from GEO accession GSE145926, read in and store as Seurat objects
nomsf <- dir("GSE145926_RAW/")[grep(".h5", dir("GSE145926_RAW/"))]
lung_obj <- list()
for(x in nomsf[1:12]){
  print(x)
  olung <- Read10X_h5(filename = paste0("GSE145926_RAW/", x))
  nom <- gsub("_filtered_feature_bc_matrix.h5", "", x)
  olung <- CreateSeuratObject(counts = olung, project = nom, min.cells = 3, min.features = 200)
  print(olung)
  lung_obj[[nom]] <- olung
}
for (i in names(lung_obj)) {
  lung_obj[[i]] <- RenameCells(lung_obj[[i]],
                               add.cell.id = i)
}
# merge all into single object and store
merged_combined <- purrr::reduce(lung_obj,merge,do.normalize = FALSE)
save(merged_combined, file="data/All_Data_combined_Raw.RData")

#### load data 
load("data/All_Data_combined_Raw.RData") # 23,742 genes x 90,696 cells
##load associated metada  
load("Source_data/metadata_scrnaseq.RData")

# Rename samples for clarity
rownames(metadata) <- metadata$geo_accession
old_sample_name <- levels(Idents(merged_combined))
sele <- NULL; for(x in old_sample_name) sele <- c(sele, grep(strsplit(x,split="_")[[1]][1], rownames(metadata)))
metadata <- metadata[sele,] 
tipus <- as.factor(metadata$patient.group.ch1); levels(tipus) <- c("healthy", "mild", "severe", "severe")
new_sample_name <- paste(unlist(strsplit(old_sample_name,"_"))[c(1:12)*2], tipus, sep="_")
levels(merged_combined@active.ident) <- new_sample_name


#####Filtering process
#####
# vcalculate percent pf mitochondria genes
##housekeeping genes
hkgenes <- read.table("Source_data/tirosh_house_keeping.txt", skip = 2)
hkgenes <- as.vector(hkgenes$V1)
merged_combined[["percent.mt"]] <- PercentageFeatureSet(merged_combined, pattern = "^(MT|mt)-")


# diagnostic plots to help with filtering
p <- ggplot(merged_combined@meta.data, aes(nCount_RNA, nFeature_RNA, color=percent.mt))
p <- p + geom_point(size=0.2)
p <- p + geom_hline(yintercept=200, color='red')
p <- p + geom_hline(yintercept=6000, color='red')


hkgenes.found <- which(toupper(rownames(merged_combined@assays$RNA@data)) %in% hkgenes)
# Add_number_of_house_keeping_genes
n.expressed.hkgenes <- Matrix::colSums(merged_combined@assays$RNA@data[hkgenes.found, ] > 0)
merged_combined <- AddMetaData(object = merged_combined, metadata = n.expressed.hkgenes, col.name = "n.exp.hkgenes")
merged_combined <- subset(merged_combined, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA > 1000 & n.exp.hkgenes > 55)



# scTransform pipeline ----------------------------------------------------

##normalize, scale and find variables in one command
merged_combined <- SCTransform(merged_combined, vars.to.regress = c("percent.mt","orig.ident"), 
                               return.only.var.gene=F,verbose = TRUE)

####PCA 
merged_combined<-RunPCA(merged_combined,verbose = T,features = )

# visualize loadings for PC1 and 2
p<-VizDimLoadings(merged_combined, dims = 1:2, reduction = "pca") # PC1: many interferon, TNFSF10 ; PC2: CD3 (T cell) markers


# UMAP
merged_combined <- RunUMAP(merged_combined, dims = 1:50) # increase threads to run faster?

# Clustering
merged_combined <- FindNeighbors(merged_combined, dims = 1:50)
merged_combined <- FindClusters(merged_combined, verbose = FALSE)
DimPlot(merged_combined, label = TRUE) + NoLegend()

### Cell type Markers

# look for cell type and relevant cell markers
# FeaturePlot(merged_combined, features = c("CD80", "CD86", "CTLA4", "IL6","MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ",
#                              "CD8A"), ncol=2)



# CD80/86 axis analysis ---------------------------------------------------

status <- tipus
sampleId <- as.factor(merged_combined@meta.data$orig.ident)
clusterAssign <- Idents(merged_combined)
sampleNames <- levels(sampleId)
Tcell_cluster <- "5" # active CD4+ T cell cluster
marker_data <- FetchData(merged_combined, c("IL6", "CD80", "CD86"))

# Number of cells per individual
TCD4_active <- IL6 <- CD80 <- CD86 <- rep(NA, length(sampleNames))
for(i in 1:length(sampleNames)){
  ix <- which(sampleId == sampleNames[i])
  # number of active T cells per sample
  clusterIds <- NULL
  clusterIds <- which(clusterAssign == Tcell_cluster)
  TCD4_active[i] <- length(intersect(ix, clusterIds))
  # number of CD80 and CD86 expressing cells per sample
  CD86[i] <- length(intersect(ix, which(marker_data[,3] > 0)))
  CD80[i] <- length(intersect(ix, which(marker_data[,2] > 0)))
  # number of IL6 producing cells per sample
  IL6[i] <- length(intersect(ix, which(marker_data[,1] > 0)))
}


# correlation analysis: test all pairwise correlations between the three elements
# of the CD80/86 axis: CD80/86+ APCs, active CD4+ T cells, and IL6 producing cells
cor.test(TCD4_active, CD86)
cor.test(TCD4_active, CD80)
cor.test(TCD4_active, IL6)
cor.test(CD86, IL6)
cor.test(CD80, IL6)


# Figure 8 plots ----------------------------------------------------------

axisData <- data.frame(sampleNames, TCD4_active, CD86, CD80, IL6, status)

fill = c("steelblue", "coral", "red")
p8A <- ggplot(axisData, aes(x = CD86 , y = TCD4_active, size=TCD4_active, fill=status)) + geom_point(shape=21) + scale_size_area(max_size = 10)
p8A <- p8A + ggtitle("") + labs(x = "#CD86 cells", y = "#Active T CD4+ cells") + scale_fill_manual(values = fill)
p8A

fill = c("steelblue", "coral", "red")
p8B <- ggplot(axisData, aes(x = CD80 , y = TCD4_active, size=TCD4_active, fill=status)) + geom_point(shape=21) + scale_size_area(max_size = 10)
p8B <- p8B + ggtitle("") + labs(x = "#CD80 cells",  y = "#Active T CD4+ cells" ) + scale_fill_manual(values = fill)
p8B

fill = c("steelblue", "coral", "red")
p8C <- ggplot(axisData, aes( x = TCD4_active, y = IL6, size=IL6, fill=status)) + geom_point(shape=21) + scale_size_area(max_size = 10)
p8C <- p8C + ggtitle("") + labs(x = "#Active T CD4+ cells", y = "#IL6+ cells") + scale_fill_manual(values = fill)
p8C

fill = c("steelblue", "coral", "red")
p8D <- ggplot(axisData, aes(x = CD86, y = IL6, size=IL6, fill=status)) + geom_point(shape=21) + scale_size_area(max_size = 10)
p8D <- p8D + ggtitle("") + labs(x = "#CD86 cells", y = "#IL6+ cells") + scale_fill_manual(values = fill)
p8D

fill = c("steelblue", "coral", "red")
p8E <- ggplot(axisData, aes(x = CD80 , y = IL6, size=IL6, fill=status)) + geom_point(shape=21) + scale_size_area(max_size = 10)
p8E <- p8E + ggtitle("") + labs(x = "#CD80 cells", y = "#IL6+ cells") + scale_fill_manual(values = fill)
p8E
