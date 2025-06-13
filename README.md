# visiumhd-helpers

R helper functions for Visium HD analysis. These functions facilitate running R-based packages used by tools like [ENACT](https://github.com/Sanofi-Public/enact-pipeline.git).

## Overview
This repository provides essential R scripts for Visium HD data analysis:
- `r_scripts/sargent.R`: Performs cell type annotation using Sargent on ENACT bin-to-cell assignment results
- `r_scripts/enact_to_seurat.R`: Converts ENACT outputs to Seurat-compatible format

## Running Sargent
### Prerequisites
Complete these ENACT pipeline steps first:
   - Segmentation
   - bin_to_geodataframes
   - bin_to_cell_assignment

Example ENACT code:
```
from enact.pipeline import ENACT

so_hd = ENACT(
    cache_dir="/home/oneai/test_cache",
    wsi_path="Visium_HD_Human_Colon_Cancer_tissue_image.btf",
    visiumhd_h5_path="binned_outputs/square_002um/filtered_feature_bc_matrix.h5",
    tissue_positions_path="binned_outputs/square_002um/spatial/tissue_positions.parquet",
    analysis_name="demo-colon",
    segmentation=True,             # ENABLE THIS STEP
    bin_to_geodataframes=True,     # ENABLE THIS STEP 
    bin_to_cell_assignment=True,   # ENABLE THIS STEP
    cell_type_annotation=False     # ** DISABLE THIS STEP **
)
```

### Installation & Setup

1. **Configure Makefile Variables**
   ```makefile
   # Update these variables in the Makefile
   ENV_DIR := /path/to/conda/environments/  # Directory to save R conda environment
   R_ENV_NAME := your_env_name             # Name for the R environment
   CONFIG_PATH := /path/to/config.yaml     # Path to ENACT config file. See `https://github.com/Sanofi-Public/enact-pipeline/blob/main/config/configs.yaml` for an example
   ```

2. **Create R Environment**

    ```make setup_r_env```

3. ***Run Sargent analysis***

    ```make run_sargent```

    This will process the ENACT results using the configuration specified in your config file.

4. ***Package Results using ENACT***
Run ENACT's cell type annotation step to:
- Package Sargent results in ENACT format
- Generate visualization .tmap file

Example ENACT code:
```
from enact.pipeline import ENACT

so_hd = ENACT(
    cache_dir="/home/oneai/test_cache",
    wsi_path="Visium_HD_Human_Colon_Cancer_tissue_image.btf",
    visiumhd_h5_path="binned_outputs/square_002um/filtered_feature_bc_matrix.h5",
    tissue_positions_path="binned_outputs/square_002um/spatial/tissue_positions.parquet",
    analysis_name="demo-colon",
    segmentation=False,             # ** DISABLE THIS STEP **
    bin_to_geodataframes=False,     # ** DISABLE THIS STEP **
    bin_to_cell_assignment=False,   # ** DISABLE THIS STEP **
    cell_type_annotation=True       # ENABLE THIS STEP
)
```