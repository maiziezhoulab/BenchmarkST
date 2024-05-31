BANKSY tutorial
============

#. Dependencies

.. code-block:: r

    .libPaths("/home/xiem6/0Virtual_Environment/R4.3_lib/")
    library(Banksy)
    library(SummarizedExperiment)
    library(SpatialExperiment)
    library(scuttle)
    library(scater)
    library(cowplot)
    library(ggplot2)
    library(Seurat)
    library(hdf5r)

#. Data loading: DLPFC

.. code-block:: r
    sample.name <- "151573"
    cluster.number <- 7
    dir.input <- file.path('/home/xiem6/0Data/DLPFC12/', sample.name)
    dir.output <- file.path('/maiziezhou_lab/manfeifei/0Projects/Benchmark_STdata/R/Banksy/output/', sample.name, '/')
    meta.input <- file.path('/home/xiem6/0Data/DLPFC12/', sample.name, 'gt')
    layer.input <- file.path('/home/xiem6/0Data/DLPFC12/', sample.name, 'gt/layered')

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
                            col.name = 'sdimx') #row' 
    sp_data <- AddMetaData(sp_data, 
                            metadata = df_meta['V4'],
                            col.name = 'sdimy') #col
    sp_data <- AddMetaData(sp_data, 
                            metadata = layer.data['layer'],
                            col.name = 'layer_guess_reordered')

    #gcm <- sp_data@assays$Spatial@counts # this is used in Seurat ‘4.3.0’
    gcm <- sp_data@assays$Spatial@layers$counts # this is used in Seurat ‘5.0.3’

    # get coordinates
    locs <- data.frame(sdimx=sp_data@meta.data$sdimx, sdimy=sp_data@meta.data$sdimy)
    row.names(locs) <- row.names(sp_data@meta.data)
    spatial_coor <- as.matrix(locs)
    rownames(spatial_coor) <- row.names(locs) 


#. Data Loading: MHypothalamus Bregma

.. code-block:: r
    dir.input <- file.path('/data/maiziezhou_lab/Datasets/ST_datasets/', sample.name)
    dir.output <- file.path('/data/maiziezhou_lab/yikang/ST_R/SpatialPCA/output/', sample.name, sheet.name)
    if(!dir.exists(file.path(dir.output))){
        dir.create(file.path(dir.output), recursive = TRUE)
    }


    filename = paste0(dir.input, '/MERFISH_Animal1_cnts.xlsx')
    cnts <- as.data.frame(read_excel(filename, sheet = sheet.name))
    row.names(cnts) <- cnts[,"...1"]
    cnts <- cnts[ -c(1) ]

    infoname = paste0(dir.input, '/MERFISH_Animal1_info.xlsx')
    xys <- as.data.frame(read_excel(infoname, sheet = sheet.name))
    row.names(xys) <- xys[,"...1"]
    xys <- xys[-c(1)]
    xy_coord <- xys[c(1,2)]

    count <- as(as.matrix(cnts), "dgCMatrix")

#. Run BANKSY

.. code-block:: r
    run_analysis <- function(respa) {
        #Initialize a SpatialExperiment object and perform basic quality control and normalization.
        se <- SpatialExperiment(assay = list(counts = gcm), spatialCoords = spatial_coor) 
        print("finish Spatial exp") 

        imgData(se) <- NULL
        assay(se, "logcounts") <- NULL
        reducedDims(se) <- NULL
        rowData(se) <- NULL
        colData(se) <- DataFrame(
            sample_id = sample.name,
            clust_annotation = factor(
                addNA(sp_data@meta.data$layer_guess_reordered),
                exclude = NULL, labels = seq(length(unique(sp_data@meta.data$layer_guess_reordered)))
            ),
            row.names = row.names(locs) 
        )
        
        rownames(se) <- row.names(sp_data) 
        #' Remove NA spots optionally
        se = se[as.numeric(se$clust_annotation) <= cluster.number, ]
        se$clust_annotation = droplevels.factor(se$clust_annotation)

        # Normalization to mean library size
        se <- computeLibraryFactors(se)
        print("finish computeLibraryFactors")
        aname = "logcounts" 
        assay(se, aname) <- normalizeCounts(se, log = TRUE)
        print("finish normalizeCounts")
        
        #' Find variable features 
        feat = VariableFeatures(FindVariableFeatures(as.Seurat(se)))
        se = se[feat, ]
        
        #Compute the neighborhood matrices for BANKSY. Setting compute_agf=TRUE computes both the weighted neighborhood mean (M
        #) and the azimuthal Gabor filter (G). The number of spatial neighbors used to compute M and G 
        #are k_geom[1]=15 and k_geom[2]=30 respectively. We run BANKSY at lambda=0 corresponding to non-spatial clustering, and lambda=0.2 corresponding to BANKSY for cell-typing.
        lambda <- c(0, 0.2)
        #k_geom <- c(15, 30) # this is default
        k_geom = c(18, 18) # this is for DLPFC
        se <- Banksy::computeBanksy(se, assay_name = aname, compute_agf = TRUE, k_geom = k_geom)
        print("finish computeBansky")

        # run PCA on the BANKSY matrix and perform clustering. Setting use_agf=TRUE uses both 
        # and to construct the BANKSY matrix.
        set.seed(1000)
        se <- Banksy::runBanksyPCA(se, use_agf = TRUE, lambda = lambda)
        print("finish PCA")
        se <- Banksy::runBanksyUMAP(se, use_agf = TRUE, lambda = lambda)
        print("finish UMAP")
        se <- Banksy::clusterBanksy(se, use_agf = TRUE, lambda = lambda, resolution = respa)
        print("finish cluster")
        #Different clustering runs can be relabeled to minimise their differences with connectClusters:
        se <- Banksy::connectClusters(se)
        print("finish connect")

        #Visualise the clustering output for non-spatial clustering (lambda=0) and BANKSY clustering (lambda=0.2).
        cnames <- colnames(colData(se))
        cnames <- cnames[grep("^clust", cnames)]
        colData(se) <- cbind(colData(se), spatialCoords(se))

        # plot_nsp <- plotColData(se,
        #     x = "sdimx", y = "sdimy",
        #     point_size = 0.6, colour_by = cnames[1]
        # )
        # plot_bank <- plotColData(se,
        #     x = "sdimx", y = "sdimy",
        #     point_size = 0.6, colour_by = cnames[2]
        # )

        #plot_grid(plot_nsp + coord_equal(), plot_bank + coord_equal(), ncol = 2)
        
        name <- paste0("clust_M1_lam0.2_k50_res", respa)
        zlabels <- colData(se)[[name]]
        # spatial domain labels, the class of zlabels is factor
        print("obatained spots no.")
        print(length(zlabels))
        ob_clusternumber <- nlevels(zlabels) # the no. of obtained clusters
        print("obtained num. of clusters")
        print(ob_clusternumber)
        gtlabels <- list(colData(se)$clust_annotation) 
        print("length of gtlables")
        print(length(gtlabels[[1]]))
        print("original num. of clusters:")
        print(length(unique(gtlabels[[1]])))
        if (ob_clusternumber == 0) {
            df_i <- data.frame()
            ari_Bansky <- 0
        } else {
            df_i <- data.frame(
            row = colData(se)$sdimx,
            col = colData(se)$sdimy,
            slot2 = zlabels,
            stringsAsFactors = FALSE
            )
            # Set row names of df_i
            rownames(df_i) <- rownames(colData(se))
            colnames(df_i)[ncol(df_i)] <- "spatial cluster"
            ari_Bansky <- mclust::adjustedRandIndex(zlabels, gtlabels[[1]]) 
        }
        
        print("ari")
        print(ari_Bansky)
        
        return(list(ari_Bansky = ari_Bansky, respa = respa, cluster.number = length(unique(gtlabels[[1]])), ob_clusternumber = ob_clusternumber, df_i=df_i, dirOut=dir.output))
    }
    # Initialize the result dataframe
    result_df <- data.frame(ari_Bansky = numeric(), res_pa = numeric(), ori_cluster_no = numeric(), ob_cluster_no = numeric())
    df_list <- list()

    # Run the main function with different resolution parameters
    respa <- seq(0.1, 1.5, by = 0.1) 

    for(i in 1:length(respa)){
    result <- run_analysis(respa[i])
    result_df <- rbind(result_df, data.frame(ari_Bansky = result$ari_Bansky, resolution_para = result$respa, ori_cluster_no = result$cluster.number, ob_cluster_no = result$ob_clusternumber))
    df_list <- c(df_list, list(result$df_i))
    }

#. Save the output

.. code-block:: r
    Write the original result dataframe to a txt file
    dir.output <- result$dirOut
    write.table(result_df, file = file.path(dir.output, "ori_ari.txt"), sep = "\t", row.names = FALSE, col.names=TRUE) 
    result_df_filtered <- result_df[result_df$ori_cluster_no == result_df$ob_cluster_no, ]
    print("filtered")
    print(result_df_filtered)
    # Check if result_df_filtered is not empty
    if (nrow(result_df_filtered) > 0) {
        # Find the row indices with ari_Bansky closest to the median value
        median_index <- which.min(abs(result_df_filtered$ari_Bansky - median(result_df_filtered$ari_Bansky)))
     
        # Extract the rows closest to the median value, just keep the first one
        median_df <- result_df_filtered[median_index[1], ]
        write.table(median_df, file = file.path(dir.output, "ari.txt"), sep = "\t", row.names = FALSE, col.names=TRUE) 
        # Print the result
        print(median_df)
        # get the corresponding row number in the original results
        para_median <- median_df$resolution_para
        print("para_median")
        print(para_median)
        row_number <- which(result_df$resolution_para == para_median)
        print("row number:")
        row_number
        df <- df_list[[row_number]]
        filename <- paste0(sample.name, "_output.csv")
        write.table(df, file = file.path(dir.output, filename), sep = "\t", qmethod = "double", col.names = NA)
    } else {
        print("Can't find a result that obtained cluster no. equal to the given cluster no.")
        lower <- result_df$ori_cluster_no - 1
        upper <- result_df$ori_cluster_no +1
        result_df_filtered <- result_df[lower <= result_df$ob_cluster_no & result_df$ob_cluster_no <= upper, ]
        print("filtered")
        print(result_df_filtered)
        if (nrow(result_df_filtered) > 0) {
            # Find the row indices with ari_Bansky closest to the median value
            median_index <- which.min(abs(result_df_filtered$ari_Bansky - median(result_df_filtered$ari_Bansky)))
            
            # Extract the rows closest to the median value, just keep the first one
            median_df <- result_df_filtered[median_index[1], ]
            write.table(median_df, file = file.path(dir.output, "ari.txt"), sep = "\t", row.names = FALSE, col.names=TRUE) 
            # Print the result
            print(median_df)
            # get the corresponding row number in the original results
            para_median <- median_df$resolution_para
            print("para_median")
            print(para_median)
            row_number <- which(result_df$resolution_para == para_median)
            print("row number:")
            row_number
            df <- df_list[[row_number]]
            filename <- paste0(sample.name, "_output.csv")
            write.table(df, file = file.path(dir.output, filename), sep = "\t", qmethod = "double", col.names = NA)
        } else {
            print("Can't find a result that obtained cluster no. in [ori_cluster_no-1, ori_cluster_no +1]")
        }
    }

