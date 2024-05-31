BASS tutorial
============

#. Dependencies

.. code-block:: r

    library(BASS)
    library(Matrix)
    library(Seurat)
    library(ggplot2)

#. Data loading: DLPFC

.. code-block:: r

    dir.input <- file.path('/data/maiziezhou_lab/Datasets/ST_datasets/DLPFC12/', sample.name)
    dir.output <- file.path('/data/maiziezhou_lab/yikang/ST_R/BASS/output/', sample.name, '/')
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
    #df_meta_matched <- df_meta[df_meta$V1 %in% row.names(sp_data@meta.data),]
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

    count <- sp_data@assays$Spatial@counts

    # get coordinates
    coord <- data.frame(row=sp_data@meta.data$row, col=sp_data@meta.data$col, annotation=sp_data@meta.data$layer_guess_reordered)
    row.names(coord) <- row.names(sp_data@meta.data)
    set.seed(0)
    # Set up BASS object
    cntm = list(count)
    xym = list(data.frame(coord[,1:2]))
    C = 20
    R = as.numeric(cluster.number)


#. Data Loading: MHypothalamus Bregma

.. code-block:: r
    
    dir.input <- file.path('/data/maiziezhou_lab/Datasets/ST_datasets/', sample.name)
    dir.output <- file.path('/data/maiziezhou_lab/yikang/ST_R/BASS/output/', sample.name, sheet.name)
    #dir.output <- file.path('/data/maiziezhou_lab/yikang/ST_R/BASS/output/', sample.name, '/')

    if(!dir.exists(file.path(dir.output))){
    dir.create(file.path(dir.output), recursive = TRUE)
    }


    filename = paste0(dir.input, '/MERFISH_Animal1_cnts.xlsx')
    cnts <- as.data.frame(read_excel(filename, sheet = sheet.name))
    row.names(cnts) <- cnts[,"...1"]
    cnts <- cnts[ -c(1) ]
    cnts <- list(cnts)

    infoname = paste0(dir.input, '/MERFISH_Animal1_info.xlsx')
    xys <- as.data.frame(read_excel(infoname, sheet = sheet.name))
    row.names(xys) <- xys[,"...1"]
    gtlabels <- list(xys$z)
    xys <- xys[-c(1)]
    xys <- xys[-c(-2:-1)]
    xys <- list(xys)

    C <- 20 # number of cell types
    R <- as.numeric(cluster.number) # number of spatial domains

#. Run BASS

.. code-block:: r

    BASS <- createBASSObject(cntm, xym, C = C, R = R, beta_method = "SW", init_method = "mclust", 
                          nsample = 1000)

    BASS <- BASS.preprocess(BASS, doLogNormalize = TRUE,
    geneSelect = "sparkx", nSE = 3000, doPCA = TRUE, 
    scaleFeature = FALSE, nPC = 20)

    # Run BASS algorithm
    BASS <- BASS.run(BASS)

    # post-process posterior samples:
    # 1.Adjust for label switching with the ECR-1 algorithm
    # 2.Summarize the posterior samples to obtain the spatial domain labels
    BASS <- BASS.postprocess(BASS)

    zlabels <- BASS@results$z # spatial domain labels


#. Calculate the ARI and save the output

.. code-block:: r

    df_i <- data.frame(slot1 = BASS@xy, slot2 = BASS@results$z)
    colnames(df_i)[ncol(df_i)] <- "spatial cluster"
    # Write the data frame to a csv file
    #write.csv(df_i, paste0("./output/data_", i, ".csv")) 
    # Loop over the numbers 1 to 4
    filename <- paste0(sample.name, "_output.csv")
    write.table(df_i, file = file.path(dir.output, filename), sep = "\t", qmethod = "double", col.names=NA)

    gtlabels <- list(sp_data@meta.data$layer_guess_reordered)
    ari_bass <- mclust::adjustedRandIndex(zlabels[[1]], gtlabels[[1]])