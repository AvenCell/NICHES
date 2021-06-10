library(SCC)
# Load Data
emb <- readRDS("/data/jyc/github_proj/dbit_seq/FFPE_paper/DBiT-seq_FFPE/Example_Data/spatial_seurat_processed_main.RDS")

#Pull out one slide only
library(Seurat)
FFPE2 <- subset(emb,idents=c("E10.5 tail"))
DefaultAssay(FFPE2) <- "RNA"
Idents(FFPE2) <- FFPE2$cell_label_dom

# Register position
FFPE2@meta.data$position <- gsub("_2","",colnames(FFPE2))

# Index position to x and y
temp <- data.frame(stringr::str_split_fixed(FFPE2$position,pattern = 'x',n = 2))
FFPE2$x <- temp$X1
FFPE2$y <- temp$X2

# Testing runs
test <- RunSCC(FFPE2,species = 'mouse',position.x = 'x',position.y = 'y',
               CellToCell = F,CellToSystem = F,SystemToCell = F,
               CellToCellSpatial =T,CellToNeighborhood = F,NeighborhoodToCell =F,meta.data.to.map = c('nCount_SCT','orig.ident')) #works

# Check with my previous implementation
library(SingleCellConnectome)
FFPE2$cell_types <- FFPE2$cell_label_dom
test1 <- SingleCellConnectome::runSCC(seu_obj = FFPE2,assay = "RNA",organizations = c("pair","pair_spatial"),metadata = c("cell_types","position"),species = "mouse",neighborhood_radius = 1,autocrine = T,n_threads = 8,downsampling = FALSE)

test2 <- SingleCellConnectome::runSCC(seu_obj = FFPE2,assay = "RNA",organizations = c("pair"),metadata = c("cell_types","position"),species = "mouse",neighborhood_radius = 1,autocrine = T,n_threads = 8,downsampling = FALSE)

test_rowsum <- rowSums(as.matrix(test[[1]]@assays$CellToCellSpatial@counts[r1,]))
test1_rowsum <- rowSums(as.matrix(test1@organizations[[2]]@assays$RNA@counts[r1,]))
sum(test_rowsum!=test1_rowsum)
