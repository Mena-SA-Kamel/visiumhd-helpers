## This vignette shows how to convert an AnnData .h5 objects from ENACT-SO to a Seurat Object
 ## Contribution by EL KHOLDI, Rayan 

 library(BPCells)
 library(Seurat)

 save_dir = "~/Documents/Bioinfo/Results"

 #We first load our AnnData object as well as the metadata file which contains both cell types and cell coordinates
 data <- open_matrix_anndata_hdf5("celltypist_results/cells_adata.h5")
 metadata <- read.csv("celltypist_results/merged_results.csv")

 #Store AnnData object as a BPCells matrix
 write_matrix_dir(mat = data,dir = save_dir)

 #Open the BPCells matrix
 h5matrix <- open_matrix_dir(dir = save_dir)

 #Create the Seurat Object and add the metadata as "cell_type" into the Seurate Object
 seurat_obj <- CreateSeuratObject(counts = h5matrix)
 seurat_obj <- AddMetaData(seurat_obj,metadata$cell_type,col.name="celltypist")

 #Creation of the Centroids and addition to the Seurat Object
 cents <- data.frame(x = metadata$cell_x,y = metadata$cell_y,cell = metadata$id)
 cents <- CreateCentroids(cents, nsides = Inf, radius = 3)
 mycoords <- list("Centroids" = cents)
 myimg <- CreateFOV(mycoords, type = c("centroids"),
                      molecules = NULL,
                      assay = "RNA",
                      key = "RNA_",
                      name = "sample")

 seurat_obj@images$sample <- myimg

 #Optional: The count matrix in the Seurat Object was created via BPCells, in that case, it's not really stored
 #in the Seurat Object but the path to the matrix is stored inside of the Seurat Object. The advantage is that
 #that the Seurat Object size is minimal. However, some downstream computations like integration won't work. We
 #can store directly the count matrix into the Seurat object with the following line of code:
 seurat_obj[["RNA"]]$counts <- as(object = seurat_obj[["RNA"]]$counts, Class = "dgCMatrix")