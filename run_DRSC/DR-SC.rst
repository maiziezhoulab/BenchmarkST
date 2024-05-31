DR-SC tutorial
============

#. Dependencies

.. code-block:: r
    
    library("DR.SC")
    library(Seurat)
    library(ggplot2)
    library(tictoc)

#. Data loading: DLPFC

.. code-block:: r

   # Loading the cell-gene expression matrix and the obs dataset for each slice.
   dir.input <- file.path('/data/maiziezhou_lab/Datasets/ST_datasets/DLPFC12/', sample.name)
   dir.output <- file.path('/data/maiziezhou_lab/yikang/ST_R/DRSC/output/', sample.name, '/')
   meta.input <- file.path('/data/maiziezhou_lab/Datasets/ST_datasets/DLPFC12/', sample.name, 'gt')
   layer.input <- file.path('/data/maiziezhou_lab/Datasets/ST_datasets/DLPFC12/', sample.name, 'gt/layered')
   #meta.input <- file.path('/data/maiziezhou_lab/yikang/ST_R/SEDR_analyses/data/DLPFC/', sample.name)
   if(!dir.exists(file.path(dir.output))){
       dir.create(file.path(dir.output), recursive = TRUE)
   }

   filename <- paste0(sample.name, "_filtered_feature_bc_matrix.h5")
   sp_data <- Load10X_Spatial(dir.input, filename = filename)

   #df_meta <- read.table(file.path(meta.input, 'metadata.tsv'))
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
                       col.name = 'annotation')

   head(sp_data)


#. Data Loading: MHypothalamus Bregma

.. code-block:: r

    dir.input <- file.path('/data/maiziezhou_lab/Datasets/ST_datasets/', sample.name)

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

#. Run the DR.SC

.. code-block:: r

    sp_data <- NormalizeData(sp_data, verbose = F)
    # choose 500 highly variable features
    seu <- FindVariableFeatures(sp_data, nfeatures = 500, verbose = F)
    ### Given K

    seu <- DR.SC(seu, K=as.numeric(cluster.number), platform = 'Visium', verbose=F)


#. Calculate the ARI

.. code-block:: r

    ## SAVE the files
    filename <- paste0(sample.name, ".csv")
    data_to_write_out <- as.data.frame(as.matrix(seu@meta.data))
    write.table(data_to_write_out, file = file.path(dir.output, filename), sep = "\t", qmethod = "double", col.names=NA)

    ## Calculate the ARI
    ari_drsc <- mclust::adjustedRandIndex(seu$spatial.drsc.cluster, seu$annotation)