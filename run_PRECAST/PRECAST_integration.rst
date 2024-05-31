PRECAST integration tutorial
============

#. Dependencies

.. code-block:: r

    library(PRECAST)
    library(Seurat)
    library(ggplot2)

#. Data loading: DLPFC

.. code-block:: r

    # Define a function to process each sample
    processSample <- function(sample.name, dir.base, dir.output.base, layer.input.base, cluster.number) {
    dir.input <- file.path(dir.base, sample.name)
    meta.input <- file.path(dir.input, 'gt')
    layer.input <- file.path(meta.input, 'layered')
    
    # Create output directory if it does not exist
    dir.output <- file.path(dir.output.base, paste0(sample.name, '_integration'))
    if(!dir.exists(dir.output)) {
        dir.create(dir.output, recursive = TRUE)
    }
    
    filename <- paste0(sample.name, "_filtered_feature_bc_matrix.h5")
    sp_data <- Load10X_Spatial(dir.input, filename = filename)
    
    df_meta <- read.table(file.path(meta.input, 'tissue_positions_list_GTs.txt'))
    split_data <- strsplit(df_meta$V1, split = ",")
    df_meta <- do.call(rbind, lapply(split_data, function(x) {
        data.frame(V1=x[1], V2=x[2], V3=as.numeric(x[3]), V4=as.numeric(x[4]), V5=x[5], V6=x[6], V7=x[7], row.names = x[1])
    }))
    
    common_cells <- colnames(sp_data[["Spatial"]]) %in% rownames(df_meta)
    sp_data <- sp_data[, common_cells]
    
    layer.data <- data.frame()
    layer.indices <- if(as.numeric(cluster.number) == 5) 3:6 else 1:6
    
    for(i in layer.indices) {
        file.name <- paste0(sample.name, "_L", i, "_barcodes.txt")
        file.path <- file.path(layer.input, file.name)
        data.temp <- read.table(file.path, header = FALSE, stringsAsFactors = FALSE)
        data.temp <- data.frame(barcode = data.temp[,1], layer = paste0("layer", i), row.names = data.temp[,1])
        layer.data <- rbind(layer.data, data.temp)
    }
    
    # For the WM file
    file.name <- paste0(sample.name, "_WM_barcodes.txt")
    file.path <- file.path(layer.input, file.name)
    data.temp <- read.table(file.path, header = FALSE, stringsAsFactors = FALSE)
    data.temp <- data.frame(barcode = data.temp[,1], layer = "WM", row.names = data.temp[,1])
    layer.data <- rbind(layer.data, data.temp)
    
    sp_data <- AddMetaData(sp_data, metadata = df_meta['V3'], col.name = 'row')
    sp_data <- AddMetaData(sp_data, metadata = df_meta['V4'], col.name = 'col')
    sp_data <- AddMetaData(sp_data, metadata = layer.data['layer'], col.name = 'layer_guess_reordered')
    
    return(sp_data)
    }

    # Define base directories
    dir.base <- '/data/maiziezhou_lab/Datasets/ST_datasets/DLPFC12/'
    dir.output.base <- '/data/maiziezhou_lab/yikang/ST_R/PRECAST/output/'
    layer.input.base <- '/gt/layered'
    sample.names <- c("sample.name1", "sample.name2") # Update these sample names accordingly
    cluster.number <- "5" # or any other number based on your requirement

    # Process each sample
    seuList <- lapply(sample.names, function(sample.name) processSample(sample.name, dir.base, dir.output.base, layer.input.base, cluster.number))


#. Data Loading: MHypothalamus Bregma

.. code-block:: r
    
    

#. Run PRECAST integration

.. code-block:: r

    set.seed(2022)
    PRECASTObj <-  CreatePRECASTObject(seuList, customGenelist=row.names(seuList[[1]]), rawData.preserve = TRUE)
    ## check the number of genes/features after filtering step
    PRECASTObj@seulist

    ## seuList is null since the default value `rawData.preserve` is FALSE.
    PRECASTObj@seuList

    ## Add adjacency matrix list for a PRECASTObj object to prepare for PRECAST model fitting.
    PRECASTObj <-  AddAdjList(PRECASTObj, platform = "Visium")

    ## Add a model setting in advance for a PRECASTObj object: verbose =TRUE helps outputing the information in the algorithm; coreNum set the how many cores are used in PRECAST. If you run PRECAST for multiple number of clusters, you can set multiple cores; otherwise, set it to 1. 
    PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal=FALSE, maxIter=30, verbose=TRUE,
                            coreNum =1)

    set.seed(2022)
    PRECASTObj <- PRECAST(PRECASTObj, K=as.numeric(cluster.number))

    ## check the fitted results: there are four list for the fitted results of each K (6:9).
    str(PRECASTObj@resList)
    ## backup the fitted results in resList
    resList <- PRECASTObj@resList
    # PRECASTObj@resList <- resList
    PRECASTObj <- SelectModel(PRECASTObj)
    ## check the best and re-organized results
    str(PRECASTObj@resList) ## The selected best K is 7

    seuInt <- IntegrateSpaData(PRECASTObj, species='Human')



#. Calculate the ARI and save the output

.. code-block:: r

    sp_data1 <- seuList[[1]]
    sp_data2 <- seuList[[2]]

    filtered_meta_data1 <- seuInt@meta.data[seuInt@meta.data$batch == 1, ]
    row.names(PRECASTObj@resList$hZ[[1]]) <- row.names(filtered_meta_data1)
    embedding1 <- PRECASTObj@resList$hZ[[1]]
    filename1 <- paste0(sample.name1, "_embeddings.csv")
    write.table(embedding1,file=file.path(dir.output, filename1), sep= "\t", qmethod = "double", col.names=NA)
    ari_precast1 <- mclust::adjustedRandIndex(filtered_meta_data1$cluster, sp_data1@meta.data$layer_guess_reordered)

    filtered_meta_data2 <- seuInt@meta.data[seuInt@meta.data$batch == 2, ]
    row.names(PRECASTObj@resList$hZ[[2]]) <- row.names(filtered_meta_data2)
    # Strip off the suffix after the `-` for both datasets.
    filtered_row_names_stripped <- sub("-.*", "", rownames(filtered_meta_data2))
    sp_data_row_names_stripped <- sub("-.*", "", row.names(sp_data2@meta.data))

    # Identify the actual names of the common cells in sp_data2
    common_cells_names <- row.names(sp_data2@meta.data)[sp_data_row_names_stripped %in% filtered_row_names_stripped]
    # Getting the indices of uncommon cells
    # Create a logical vector indicating whether each element of sp_data_row_names_stripped is in filtered_row_names_stripped
    common_indices <- sp_data_row_names_stripped %in% filtered_row_names_stripped

    # Negate the common_indices vector to get a logical vector for uncommon indices
    uncommon_indices <- !common_indices

    # Now, extract the row names from sp_data2@meta.data that are not present in common_cells_names
    uncommon_row_names <- row.names(sp_data2@meta.data)[uncommon_indices]

    # Printing the uncommon row names
    print(uncommon_row_names)
    # Subset sp_data2 using the names of the common cells
    sp_data2 <- sp_data2[, common_cells_names]
    #common_cells2 <- colnames(sp_data2[["Spatial"]]) %in% rownames(filtered_meta_data2)
    # Subset sp_data to keep only these cells
    #filtered_meta_data2 <- filtered_meta_data2[, common_cells2]
    embedding2 <- PRECASTObj@resList$hZ[[2]]
    filename2 <- paste0(sample.name2,"_embeddings.csv")
    write.table(embedding2,file=file.path(dir.output, filename2), sep= "\t", qmethod = "double", col.names=NA)
    ari_precast2 <- mclust::adjustedRandIndex(filtered_meta_data2$cluster, sp_data2@meta.data$layer_guess_reordered)


    #cluster_df1 <- seuInt@meta.data[seuInt@meta.data$batch == 1, "cluster", drop = FALSE]
    cluster_df1 <- cbind(sp_data1@meta.data, seuInt@meta.data[seuInt@meta.data$batch == 1, "cluster", drop = FALSE])
    filename3 <- paste0(sample.name1, "_cluster.csv")
    write.table(cluster_df1,file=file.path(dir.output, filename3), sep= "\t", qmethod = "double", col.names=NA)

    #cluster_df2 <- seuInt@meta.data[seuInt@meta.data$batch == 2, "cluster", drop = FALSE]
    cluster_df2 <- cbind(sp_data2@meta.data, seuInt@meta.data[seuInt@meta.data$batch == 2, "cluster", drop = FALSE])
    filename4 <- paste0(sample.name2, "_cluster.csv")
    write.table(cluster_df2,file=file.path(dir.output, filename4), sep= "\t", qmethod = "double", col.names=NA)