# SingleCellConnectome (SCC)

### Download and Install
   
    library(devtools)
   
    install_github("KevinBastianYang/SingleCellConnectome",auth_token = "3c8b47b19256505ef76d1046b25e1dc817d87ac8",dependencies =TRUE)

### Run example: as in test/test_new.R

    library(SingleCellConnectome)
   ##### load the (spatial) scrna-seq data: change the path to yours
    spatial_data <- readRDS("./data/spatial_test_data.RDS")

   ##### Run SCC
    scc_obj <- runSCC(seu_obj = spatial_data,organizations=c("pair","pair_spatial","niche_spatial"),
                               assay='alra',metadata = c("seurat_clusters","position","cluster"),neighborhood_radius = 1,
                               species = "mouse",downsampling = FALSE,sample_size_n = 10,autocrine =TRUE,n_threads = 8)
   ##### Downstream functions                 
    scc_obj <- FindVariableFeatures(scc_obj, organization = "pair_spatial",selection.method = "vst")
    scc_obj <- ScaleData(scc_obj, organization = "pair_spatial")
    scc_obj <- RunPCA(scc_obj, organization = "pair_spatial",features = VariableFeatures(object = scc_obj,organization = "pair_spatial"),seed.use = 34)
    scc_obj <- RunTSNE(object=scc_obj,organization = "pair_spatial",reduction = "pca",dims= 1:20,seed.use=as.integer(34))
    DimPlot(object=scc_obj,organization = "pair_spatial",reduction = "tsne",group.by="cp_types")

    Seurat::Idents(scc_obj@organizations$pair_spatial) <- scc_obj@organizations$pair_spatial@meta.data$cp_types

    FeaturePlot(object=scc_obj,organization = "pair_spatial",features=c("Adam10-Epha3","Adam12-Itga9","Adam9-Itga6"))

    DotPlot(object=scc_obj,organization = "pair_spatial",features=c("Adam10-Epha3","Adam12-Itga9","Adam9-Itga6"))
