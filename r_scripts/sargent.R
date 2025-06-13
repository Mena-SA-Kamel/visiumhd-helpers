# Script to run Sargent on the output bin to cell assignment results

# Load required packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(sargent)
library(Matrix)
library(yaml)
library(dplyr)
library(readr)

# Define the list of files
# configs_file <- "/home/oneai/visiumhd-helpers/sample_configs.yaml" # Replace with path to your configurations file

# Read configs_file from environment variable
configs_file <- Sys.getenv("CONFIG_PATH")

# Add a check to ensure the environment variable is set
if (configs_file == "") {
  stop("Environment variable CONFIG_PATH is not set. Please set it before running this script.")
}


# Load the YAML file
configs <- yaml.load_file(configs_file)
# Convert the loaded data into the desired format
gene.sets <- configs$cell_markers

# Extract the directory variable from the YAML content
analysis_name <- configs$analysis_name
bin_to_cell_method <- configs$params$bin_to_cell_method
cache_dir <- configs$cache_dir

# Construct the input directory path using the variable
input_dir <- file.path(cache_dir, analysis_name, "chunks", bin_to_cell_method, "bin_to_cell_assign")
output_dir <- file.path(cache_dir, analysis_name, "chunks", bin_to_cell_method, "sargent_results")
cell_ix_lookup_dir <- file.path(cache_dir, analysis_name, "chunks", bin_to_cell_method, "cell_ix_lookup")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
# List all files in the input directory
files <- list.files(input_dir, pattern="*.csv", full.names=TRUE)

# Loop through each file
for (file in files) {
    # Skip the file if the output file already exists
    # Define the output file path
    output_file <- file.path(output_dir, basename(file))
    
    if (file.exists(output_file)) {
        cat("Skipping:", file, "(output file already exists)\n")
        next
    }
    cat("Started:", file, "\n")
    # Read the data
    countsData <- read.csv(file, header=TRUE)
    dataMatrix <- as.matrix(countsData)
    sparse_Matrix <- as(dataMatrix, "sparseMatrix")
    gex <- t(sparse_Matrix)
    colnames(gex) <- unlist(gex[1, ])
    gex <- gex[-1, ]
    tryCatch({
        adjacent.mtx <- attr(CreateSeuratObject(counts=gex) %>%
                            NormalizeData(., normalization.method="LogNormalize", 
                                        scale.factor=1e6, verbose=FALSE) %>%
                            FindVariableFeatures(., selection.method="vst", 
                                                nfeatures=2000, verbose=FALSE) %>%
                            ScaleData(., do.scale=TRUE, do.center=TRUE, 
                                    verbose=FALSE) %>%
                            RunPCA(., features=VariableFeatures(.), 
                                    verbose=FALSE) %>%
                            FindNeighbors(., reduction="pca", dims=1:30, 
                                        k.param=20, verbose=FALSE), 
                        which="graphs")[["RNA_nn"]]
        srgnt <- sargentAnnotation(gex=gex,
                            gene.sets=gene.sets,
                            adjacent.mtx=adjacent.mtx)
        
        predicted <- fetchAssignment(srgnt)
        write.csv(predicted, file=output_file, row.names=TRUE)
        cat("Processed:", file, "\n")
        
    }, error = function(e) {
    message("An error occurred: ", e$message)
    return(NULL)
    })
    }


# List to store results from all patches
sargent_results_list <- list()
chunks <- list.files(output_dir)

# Iterate over each chunk
for (chunk_name in chunks) {
    # Skip specific files
    if (chunk_name %in% c("merged_results.csv", ".ipynb_checkpoints", "eval")) {
        next
        }
    # Read the CSV files
    cell_labels <- read_csv(file.path(output_dir, chunk_name), show_col_types = FALSE)
    index_lookup <- read_csv(file.path(cell_ix_lookup_dir, chunk_name), show_col_types = FALSE)

    # Combine the data frames
    sargent_result_chunk <- bind_cols(index_lookup, cell_labels["x"])
    # sargent_result_chunk <- select(sargent_result_chunk, -`Unnamed: 0`) /

    # Add the result chunk to the list
    sargent_results_list[[length(sargent_results_list) + 1]] <- sargent_result_chunk
}
# Combine all result chunks into a single data frame
sargent_results_df <- bind_rows(sargent_results_list)

print(sargent_results_df)
# Rename the column
sargent_results_df <- rename(sargent_results_df, cell_type = x)

# Write the final result to a CSV file
write_csv(sargent_results_df, file.path(output_dir, "merged_results.csv"))