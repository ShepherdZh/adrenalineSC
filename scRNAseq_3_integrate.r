library(Seurat)

organoid1_mtx=Read10X('/pathto/organoid1/filtered_feature_bc_matrix')
organoid2_mtx=Read10X('/pathto/organoid2/filtered_feature_bc_matrix')
organoid3_mtx=Read10X('/pathto/organoid3/filtered_feature_bc_matrix')

organoid1=CreateSeuratObject(organoid1_mtx)
organoid2=CreateSeuratObject(organoid2_mtx)
organoid3=CreateSeuratObject(organoid3_mtx)

organoid1$orig.ident='organoid1'
organoid2$orig.ident='organoid2'
organoid3$orig.ident='organoid3'

mtx=read.csv('pathto/mtx.csv',row.names = 1)
ann=read.csv('pathto/obs.csv')
genes=read.csv('pathto/genes.csv')
cells=read.csv('pathto/cells.csv')

colnames(mtx) <- genes$X0
rownames(mtx) <- cells$X0

mtx <- t(mtx)

tissue <- CreateSeuratObject(mtx)

tissue$celltype <- ann$PAGODA_hc
tissue$sample <- ann$samples

split_seurat <- list(organoid1=organoid1,organoid2=organoid2,organoid3=organoid3,tissue=tissue)

split_seurat <- lapply(X = split_seurat, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

integ_features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000) 

integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, normalization.method = "LogNormalize", anchor.features = integ_features, dims = 1:40, reduction = "cca")

seurat_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "LogNormalize", dims = 1:40)

seurat_integrated <- ScaleData(seurat_integrated, verbose = FALSE)
seurat_integrated <- RunPCA(object = seurat_integrated)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)
seurat_integrated <- RunTSNE(seurat_integrated, dims = 1:40)
seurat_integrated <- FindNeighbors(object = seurat_integrated, dims = 1:40)
seurat_integrated <- FindClusters(object = seurat_integrated, resolution = c(0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5))

DefaultAssay(x) <- 'RNA'
seurat_integrated$mitoratio <- PercentageFeatureSet(seurat_integrated,pattern = "^MT-")
seurat_integrated$mitoratio=(seurat_integrated$mitoratio)/100

saveRDS(seurat_integrated,'adrenaline_tis_org.rds')
