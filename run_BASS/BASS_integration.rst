BASS integration tutorial
============

#. Dependencies

.. code-block:: r

    library(BASS)
    library(Matrix)
    library(Seurat)
    library(ggplot2)

#. Data loading: DLPFC

.. code-block:: r

    # Define a function to load spatial data and process metadata
    load_and_process_sample <- function(sample_name, dir_base, dir_meta, dir_layer) {
    # Construct file paths
    dir_input <- file.path(dir_base, sample_name)
    meta_input <- file.path(dir_meta, sample_name, 'gt')
    layer_input <- file.path(dir_meta, sample_name, 'gt/layered')

    # Load spatial data
    filename <- paste0(sample_name, "_filtered_feature_bc_matrix.h5")
    sp_data <- Load10X_Spatial(dir_input, filename = filename)

    # Process and add metadata
    df_meta <- read.table(file.path(meta_input, 'tissue_positions_list_GTs.txt'))
    split_data <- strsplit(df_meta$V1, split = ",")
    df_meta <- do.call(rbind, lapply(split_data, function(x) {
    data.frame(V1=x[1], V2=x[2], V3=as.numeric(x[3]), V4=as.numeric(x[4]), V5=x[5], V6=x[6], V7=x[7])
    }))
    row.names(df_meta) <- df_meta$V1

    common_cells <- colnames(sp_data[["Spatial"]]) %in% rownames(df_meta)
    sp_data <- sp_data[, common_cells]

    # Initialize dataframe for layer data
    layer_data <- data.frame()

    # Determine layer range based on cluster number
    layer_range <- if(as.numeric(cluster.number) == 5) 3:6 else 1:6

    # Process layer data
    for(i in layer_range) {
    file_name <- paste0(sample_name, "_L", i, "_barcodes.txt")
    data_temp <- read.table(file.path(layer_input, file_name), header = FALSE, stringsAsFactors = FALSE)
    data_temp <- data.frame(barcode = data_temp[,1], layer = paste0("layer", i), row.names = data_temp[,1])
    layer_data <- rbind(layer_data, data_temp)
    }

    # Process WM layer
    file_name <- paste0(sample_name, "_WM_barcodes.txt")
    data_temp <- read.table(file.path(layer_input, file_name), header = FALSE, stringsAsFactors = FALSE)
    data_temp <- data.frame(barcode = data_temp[,1], layer = "WM", row.names = data_temp[,1])
    layer_data <- rbind(layer_data, data_temp)

    # Add metadata to spatial data
    sp_data <- AddMetaData(sp_data, metadata = df_meta['V3'], col.name = 'row')
    sp_data <- AddMetaData(sp_data, metadata = df_meta['V4'], col.name = 'col')
    sp_data <- AddMetaData(sp_data, metadata = layer_data['layer'], col.name = 'layer_guess_reordered')

    return(sp_data)
    }

    # Base directory paths
    dir_base <- '/data/maiziezhou_lab/Datasets/ST_datasets/DLPFC12/'
    dir_meta_base <- '/data/maiziezhou_lab/Datasets/ST_datasets/DLPFC12/'
    dir_layer_base <- '/data/maiziezhou_lab/Datasets/ST_datasets/DLPFC12/'
    dir_output_base <- '/data/maiziezhou_lab/yikang/ST_R/BASS/output/'

    # Sample names
    sample_names <- c("sample.name1", "sample.name2")

    # Loop through samples
    for(sample_name in sample_names) {
    dir_output <- file.path(dir_output_base, paste0(sample_name, '_integration'))
    if(!dir.exists(dir_output)) {
    dir.create(dir_output, recursive = TRUE)
    }

    # Call function to load and process each sample
    sp_data <- load_and_process_sample(sample_name, dir_base, dir_meta_base, dir_layer_base)

    }

    counts_list <- list()
    coords_list <- list()

    # Loop through each spatial data object to extract counts and coordinates
    for (i in 1:length(sample_names)) {
    # Assuming sp_data_list is a list of spatial data objects returned by load_and_process_sample
    sp_data <- sp_data_list[[i]]

    # Extract counts
    counts_list[[i]] <- sp_data@assays$Spatial@counts

    # Extract coordinates and format them
    coords <- data.frame(row = sp_data@meta.data$row, col = sp_data@meta.data$col)
    row.names(coords) <- row.names(sp_data@meta.data)
    coords_list[[i]] <- coords
    }

    # Assuming cluster.number is defined elsewhere in the script
    R <- as.numeric(cluster.number)

    # Set a fixed seed for reproducibility
    set.seed(0)

    # Parameters for BASS
    C <- 20  # Number of clusters, for example


#. Data Loading: MHypothalamus Bregma

.. code-block:: r
    
    

#. Run BASS integration

.. code-block:: r

    BBASS <- createBASSObject(cntm, xym, C = C, R = R, beta_method = "SW", init_method = "mclust", 
                      nsample = 1000)

    BASS <- BASS.preprocess(BASS, doLogNormalize = TRUE, geneSelect = "sparkx", nSE = 3000, doPCA = TRUE, scaleFeature = FALSE, nPC = 20)

    # Run BASS algorithm
    BASS <- BASS.run(BASS)

    # post-process posterior samples:
    # 1.Adjust for label switching with the ECR-1 algorithm
    # 2.Summarize the posterior samples to obtain the spatial domain labels
    BASS <- BASS.postprocess(BASS)


#. Calculate the ARI and save the output

.. code-block:: r

for (i in 1:length(sample_names)) {
    sample_name <- sample_names[i]

    # Extract spatial domain labels and coordinates from BASS results
    zlabels <- BASS@results$z[[i]]  # Spatial domain labels for current sample
    coords <- BASS@xy[[i]]  # Coordinates for current sample

    # Create a data frame for output
    df_output <- data.frame(coords, spatial_cluster = zlabels)
    colnames(df_output)[1:2] <- c("row", "col")  # Assuming coords is a dataframe with row and col

    # Save the output data frame to a CSV file
    filename <- paste0(sample_name, "_output.csv")
    write.table(df_output, file = file.path(dir_output_base, paste0(sample_name, '_integration'), filename), sep = "\t", quote = FALSE, row.names = FALSE)

    # Calculate ARI using ground truth labels
    # Assuming ground truth labels (layer_guess_reordered) are stored within the meta.data of sp_data objects
    gtlabels <- sp_data_list[[i]]@meta.data$layer_guess_reordered
    ari_bass <- mclust::adjustedRandIndex(zlabels, gtlabels)

    # Print ARI for each sample
  cat("ARI for", sample_name, ":", ari_bass, "\n")
}