PRECAST tutorial
============

#. Dependencies

.. code-block:: r

    library(PRECAST)
    library(Seurat)
    library(ggplot2)

#. Data loading: DLPFC

.. code-block:: r

    dir.input <- file.path('/data/maiziezhou_lab/Datasets/ST_datasets/DLPFC12/', sample.name)
    dir.output <- file.path('/data/maiziezhou_lab/yikang/ST_R/PRECAST/output/', sample.name, '/')
    meta.input <- file.path('/data/maiziezhou_lab/Datasets/ST_datasets/DLPFC12/', sample.name, 'gt')

    layer.input <- file.path('/data/maiziezhou_lab/Datasets/ST_datasets/DLPFC12/', sample.name, 'gt/layered')
    
    if(!dir.exists(file.path(dir.output))){
        dir.create(file.path(dir.output), recursive = TRUE)
    }

    filename <- paste0(sample.name, "_filtered_feature_bc_matrix.h5")
    sp_data <- Load10X_Spatial(dir.input, filename = filename)


    df_meta <- read.table(file.path(meta.input, 'tissue_positions_list_GTs.txt'))


    original_row_names <- row.names(df_meta) 
    split_data <- strsplit(df_meta$V1, split = ",")
    df_meta <- do.call(rbind, lapply(split_data, function(x) {
    data.frame(V1=x[1], V2=x[2], V3=x[3], V4=x[4], V5=x[5], V6=x[6], V7=x[7])
    }))
    row.names(df_meta) <- df_meta$V1
    df_meta$V3 <- as.numeric(df_meta$V3)
    df_meta$V4 <- as.numeric(df_meta$V4)

    # Set the row names of df_meta_matched to be V1
    # Identify the cells that are in both sp_data and df_meta
    common_cells <- colnames(sp_data[["Spatial"]]) %in% rownames(df_meta)

    # Subset sp_data to keep only these cells
    sp_data <- sp_data[, common_cells]

    # Initialize an empty dataframe to hold the final results
    layer.data <- data.frame()

    if(as.numeric(cluster.number) == 5) {
    for(i in 3:6){
        file.name <- paste0(sample.name, "_L", i, "_barcodes.txt")
        file.path <- file.path(layer.input, file.name)

        data.temp <- read.table(file.path, header = FALSE, stringsAsFactors = FALSE) # assuming the file has no header
        data.temp <- data.frame(barcode = data.temp[,1], layer = paste0("layer", i), row.names = data.temp[,1])

        # Append to the final dataframe
        layer.data <- rbind(layer.data, data.temp)
    }
    } else {
    for(i in 1:6){
        file.name <- paste0(sample.name, "_L", i, "_barcodes.txt")
        file.path <- file.path(layer.input, file.name)

        data.temp <- read.table(file.path, header = FALSE, stringsAsFactors = FALSE) # assuming the file has no header
        data.temp <- data.frame(barcode = data.temp[,1], layer = paste0("layer", i), row.names = data.temp[,1])

        # Append to the final dataframe
        layer.data <- rbind(layer.data, data.temp)
    }
    }


    # For the WM file
    file.name <- paste0(sample.name, "_WM_barcodes.txt")
    file.path <- file.path(layer.input, file.name)

    data.temp <- read.table(file.path, header = FALSE, stringsAsFactors = FALSE) # assuming the file has no header
    data.temp <- data.frame(barcode = data.temp[,1], layer = "WM", row.names = data.temp[,1])

    # Append to the final dataframe
    layer.data <- rbind(layer.data, data.temp)



    sp_data <- AddMetaData(sp_data, 
                        metadata = df_meta['V3'],
                        col.name = 'row')
    sp_data <- AddMetaData(sp_data, 
                        metadata = df_meta['V4'],
                        col.name = 'col')
    sp_data <- AddMetaData(sp_data, 
                        metadata = layer.data['layer'],
                        col.name = 'layer_guess_reordered')

    head(sp_data@meta.data)


#. Data Loading: MHypothalamus Bregma

.. code-block:: r
    
    dir.input <- file.path('/data/maiziezhou_lab/Datasets/ST_datasets/', sample.name)
    dir.output <- file.path('/data/maiziezhou_lab/yikang/ST_R/PRECAST/output/', sample.name, sheet.name)

    if(!dir.exists(file.path(dir.output))){
        dir.create(file.path(dir.output), recursive = TRUE)
    }


    filename = paste0(dir.input, '/MERFISH_Animal1_cnts.xlsx')
    cnts <- as.data.frame(read_excel(filename, sheet = sheet.name))
    row.names(cnts) <- cnts[,"...1"]
    cnts <- cnts[ -c(1) ]
    #cnts <- list(cnts)

    infoname = paste0(dir.input, '/MERFISH_Animal1_info.xlsx')
    xys <- as.data.frame(read_excel(infoname, sheet = sheet.name))
    row.names(xys) <- xys[,"...1"]
    xys <- xys[-c(1)]

    sp_data <- CreateSeuratObject(counts = cnts, project = "43F", min.cells = 3, names.delim = "-", names.field = 2)

    sp_data <- AddMetaData(sp_data, 
                    metadata = xys$x,
                    col.name = 'row')
    sp_data <- AddMetaData(sp_data, 
                    metadata = xys$y,
                    col.name = 'col')
    sp_data <- AddMetaData(sp_data, 
                    metadata = xys$z,
                    col.name = 'layer_guess_reordered')

    sp_data$orig.ident <- 1
    Idents(sp_data) <- row.names(sp_data@meta.data)

#. Run PRECAST clustering

.. code-block:: r

    set.seed(2023)
    preobj <- CreatePRECASTObject(seuList = list(sp_data), selectGenesMethod = "HVGs", gene.number = 2000)  # 

    preobj@seulist
    PRECASTObj <- AddAdjList(preobj, platform = "Visium")
    PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = FALSE, coreNum = 1, maxIter = 30, verbose = TRUE)

    PRECASTObj <- PRECAST(PRECASTObj, K = as.numeric(cluster.number))


#. Calculate the ARI and save the output

.. code-block:: r

    resList <- PRECASTObj@resList
    PRECASTObj <- SelectModel(PRECASTObj)
    ari_precast <- mclust::adjustedRandIndex(PRECASTObj@resList$cluster[[1]], PRECASTObj@seulist[[1]]$layer_guess_reordered)

    seuInt <- PRECASTObj@seulist[[1]]
    seuInt@meta.data$cluster <- factor(unlist(PRECASTObj@resList$cluster))
    seuInt@meta.data$batch <- 1
    seuInt <- Add_embed(PRECASTObj@resList$hZ[[1]], seuInt, embed_name = "PRECAST")
    row.names(PRECASTObj@resList$hZ[[1]]) <- row.names(seuInt@meta.data)
    embedding <- PRECASTObj@resList$hZ[[1]]
    filename <- paste0(sample.name, "_embeddings.csv")
    write.table(embedding,file=file.path(dir.output, filename), sep= "\t", qmethod = "double", col.names=NA)

    posList <- lapply(PRECASTObj@seulist, function(x) cbind(x$row, x$col))
    seuInt <- Add_embed(posList[[1]], seuInt, embed_name = "position")
    Idents(seuInt) <- factor(seuInt@meta.data$cluster)

    filename <- paste0(sample.name, "_output.csv")
    data_to_write_out <- as.data.frame(as.matrix(seuInt@meta.data))
    write.table(data_to_write_out, file = file.path(dir.output, filename), sep = "\t", qmethod = "double", col.names=NA)