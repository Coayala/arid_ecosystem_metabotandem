############################
# PRE-PROCESSING FUNCTIONS #
############################

# -------------------------------------------------------------------------
load_spectra_data <- function(datadir, metadata_file, format = 'mzML'){
  ##### Load spectra data ####
  
  # Open metadata
  metadata <- read_csv(metadata_file)
  
  # Get list of mass spectra files
  ms_files <- list.files(data_dir, full.names = TRUE, pattern = paste0('*.', format))
  
  # Read data as an `OnDiskMSnExp` object from xcms
  data <- readMSData(ms_files, 
                     pdata = new('NAnnotatedDataFrame', metadata),
                     mode = 'onDisk', 
                     verbose = TRUE)
  res <- list(spectra_data = data,
              metadata = metadata)
  
  return(res)
}

# -------------------------------------------------------------------------
centroid_check <- function(spectra_data, transform = TRUE){
  #### Check if data is centroided and transform if needed ####
  
  is.centroided <- unique(fData(data$spectra_data)$centroided)
  if(is.centroided){
    print('Data is centroided')
  } else{
    print('Data is not centroided')
    if(transform){
      print('Transforming data')
      data_cent <- spectra_data %>% 
        smooth(method = "SavitzkyGolay") %>% 
        pickPeaks()
    }
  }
  
  return(data_cent)
}

# -------------------------------------------------------------------------
calculate_tic <- function(data, show.plot = TRUE, save.figure = FALSE, color_by = 'AUTO'){
  #### Calculate, plot and save total ion chromatogram
  
  tic <- chromatogram(data$spectra_centroid, aggregationFun="sum")
  
  if(color_by == 'AUTO'){
    color_vector <- ggpubr::get_palette(palette = 'd3', length(data$metadata$SampleID))
    names(color_vector) <- data$metadata$SampleID
  } else{
    treatments <- pull(data$metadata, color_by)
    colors <- ggpubr::get_palette(palette = 'd3', length(unique(treatments)))
    names(colors) <- unique(treatments)
    
    color_vector <- c()
    for(i in 1:length(pull(data$metadata, color_by))){
      color_vector[i] <- colors[treatments[i]]
      names(color_vector)[i] <- treatments[i]
    }
    
    
  }
  
  tic_plot <- plot(tic, col = color_vector, ylab = "Intensity", xlab = "Retention Time (sec)", main = 'TIC')
  
  if(show.plot){
    tic_plot
    legend("topright", legend = unique(names(color_vector)), col = unique(color_vector), lty=1)
  }
  
  if(save.figure == FALSE){
    print('TIC plot will not be saved')
  } else{
    png(save.figure, width = 3000, height = 2000, res = 300)
    plot(tic, col = color_vector, ylab = "Intensity", xlab = "Retention Time (sec)", main = 'TIC')
    legend("topright", legend = unique(names(color_vector)), col = unique(color_vector), lty=1)
    dev.off()
  }
  
  res <- list(tic = tic,
              tic_plot = tic_plot)
  
  return(res)
}

# -------------------------------------------------------------------------
calculate_bpc <- function(data, show.plot = TRUE, save.figure = FALSE, color_by = 'AUTO'){
  #### Calculate, plot and save total ion chromatogram
  
  bpc <- chromatogram(data$spectra_centroid, aggregationFun="max")
  
  if(color_by == 'AUTO'){
    color_vector <- ggpubr::get_palette(palette = 'd3', length(data$metadata$SampleID))
    names(color_vector) <- data$metadata$SampleID
  } else{
    treatments <- pull(data$metadata, color_by)
    colors <- ggpubr::get_palette(palette = 'd3', length(unique(treatments)))
    names(colors) <- unique(treatments)
    
    color_vector <- c()
    for(i in 1:length(pull(data$metadata, color_by))){
      color_vector[i] <- colors[treatments[i]]
      names(color_vector)[i] <- treatments[i]
    }
    
    
  }
  
  bpc_plot <- plot(bpc, col = color_vector, ylab = "Intensity", xlab = "Retention Time (sec)", main = 'BPC')
  
  if(show.plot){
    bpc_plot
    legend("topright", legend = unique(names(color_vector)), col = unique(color_vector), lty=1)
  }
  
  if(save.figure == FALSE){
    print('BPC plot will not be saved')
  } else{
    png(save.figure, width = 3000, height = 2000, res = 300)
    plot(bpc, col = color_vector, ylab = "Intensity", xlab = "Retention Time (sec)", main = 'BPC')
    legend("topright", legend = unique(names(color_vector)), col = unique(color_vector), lty=1)
    dev.off()
  }
  
  res <- list(bpc = bpc,
              bpc_plot = bpc_plot)
  
  return(res)
}

# -------------------------------------------------------------------------
test_peak_picking <- function(data_centroid, mz.range, rt.range, p.width, snt, noise){
  # Test peak picking parameters
  
  cwp <- CentWaveParam(
    peakwidth = p.width,
    snthresh = snt,
    noise = noise,
    prefilter = c(1, 100)
  )
  
  data_centroid %>% 
    filterRt(rt = rt.range) %>% 
    filterMz(mz = mz.range) %>% 
    chromatogram(., aggregationFun="max") %>% 
    findChromPeaks(., param = cwp) %>% 
    plot(col = "indianred2", 
         ylab="Intensity", xlab="Retention Time (sec)",
         font.lab=1, cex.lab=1, cex.axis=1, font.main=1, cex.main=1)
  
  
}

# -------------------------------------------------------------------------
apply_peak_picking <- function(data_centroid, p.width, snt, noise){
  # Test peak picking parameters
  
  cwp <- CentWaveParam(
    peakwidth = p.width,
    snthresh = snt,
    noise = noise,
    prefilter = c(1, 100)
  )
  
  spectra_data <- findChromPeaks(data$spectra_centroid, param = cwp)
  
  return(spectra_data)
}

# -------------------------------------------------------------------------
extract_number_peaks <- function(data_centroid, save.table = FALSE, save.figure = FALSE){
  chrom_peaks_df <- as.data.frame(chromPeaks(data_centroid)) %>% 
    count(sample) %>% 
    rename('sample_index' = 'sample', 'total_peaks_detected' = 'n') %>% 
    left_join(rowid_to_column(pData(data_centroid)), by = c('sample_index' = 'rowid')) %>% 
    select(-FileName)
  
  par(mar=c(5,6,1,1))
  plotChromPeakImage(data_centroid, 
                     binSize=100, 
                     xlab="Retention time (sec)",
                     cex.sub = 0.5,
                     yaxt = "n",
                     main = 'Numbers of peaks detected')
  axis(2, 
       at = seq(0,1, length.out = length(pData(data_centroid)$SampleID)), 
       labels = pData(data_centroid)$SampleID, 
       cex.axis = 1,
       las = 2)
  
  if(save.table == FALSE){
    print('Table not saved')
  } else {
    write_csv(chrom_peaks_df, save.table)
  }
  
  if(save.figure == FALSE){
    print('Figure not saved')
  } else {
    png(save.figure, width = 3000, height = 2000, res = 300)
    colfunc <- colorRampPalette(c("lightblue", "blue4"))
    par(mar=c(5,6,1,1))
    plotChromPeakImage(data_centroid, 
                       binSize=100, 
                       xlab="Retention time (sec)",
                       cex.sub = 0.5,
                       yaxt = "n",
                       main = 'Numbers of peaks detected')
    axis(2, 
         at = seq(0,1, length.out = length(pData(data_centroid)$SampleID)), 
         labels = pData(data_centroid)$SampleID, 
         cex.axis = 1,
         las = 2)
    dev.off()
  }
  
  return(chrom_peaks_df)
}


# -------------------------------------------------------------------------
apply_alignment <- function(data, minFraction, minSamples, binSize, group_by){
  # Defining peak density parameters
  
  sample_groups <- pull(data$metadata, group_by)
  pdp <- PeakDensityParam(sampleGroups = sample_groups, 
                          bw = 30, 
                          minFraction = minFraction, 
                          minSamples = minSamples,
                          binSize = binSize)
  
  # Defining peak grouping parameters
  pgp <- PeakGroupsParam(minFraction = minFraction)
  
  ## a - Group peaks
  data_grouped <- groupChromPeaks(data$spectra_centroid, param = pdp)
  ## b - alignment
  data_aligned <- adjustRtime(data_grouped, param = pgp)
  
  treatments <- pull(data$metadata, group_by)
  colors <- ggpubr::get_palette(palette = 'd3', length(unique(treatments)))
  names(colors) <- unique(treatments)
  
  color_vector <- c()
  for(i in 1:length(pull(data$metadata, group_by))){
    color_vector[i] <- colors[treatments[i]]
    names(color_vector)[i] <- treatments[i]
  }
  
  plotAdjustedRtime(data_aligned,
                    col = color_vector, 
                    xlab="Retention Time (sec)", 
                    font.lab=2, cex.lab=2, cex.axis=2, 
                    font.main=2, cex.main=2, lwd=2)
  legend("topright", legend = unique(names(color_vector)), col = unique(color_vector), lty=1)
  
  return(data_aligned)
}

# -------------------------------------------------------------------------
apply_correspondence <- function(data, minFraction, minSamples, binSize, group_by){
  # Defining peak density parameters
  
  sample_groups <- pull(data$metadata, group_by)
  pdp <- PeakDensityParam(sampleGroups = sample_groups, 
                          bw = 30, 
                          minFraction = minFraction, 
                          minSamples = minSamples,
                          binSize = binSize)
  
  ## a - Group peaks
  data_grouped <- groupChromPeaks(data$spectra_centroid, param = pdp)
  
  return(data_grouped)
  
}

# -------------------------------------------------------------------------
apply_gap_filling <- function(data){
  
  ## determine the number of missing values
  number_na_i = sum(is.na(featureValues(data$spectra_centroid)))
  
  ## a - define parameter
  fpp <- FillChromPeaksParam(ppm = 2, expandMz = 0.25)
  
  ## b - fill in
  data_gap_filled <- fillChromPeaks(data$spectra_centroid, param=fpp)
  
  ## remaining number of na values
  number_na_f = sum(is.na(featureValues(data_gap_filled)))
  print(paste('The number of gap filled peaks was', number_na_i - number_na_f))
  
  return(data_gap_filled)
}

# -------------------------------------------------------------------------
extract_features <- function(data, save.table = FALSE){
  
  ## extract feature values after filling in
  feature_abundance_matrix <- as.data.frame(featureValues(data$spectra_centroid, value="into", method="maxint")) %>% 
    rownames_to_column(var = 'FeatureID')
  ## replace NA with zero
  feature_abundance_matrix[is.na(feature_abundance_matrix)] <- 0
  ## replace file name with sample name
  colnames(feature_abundance_matrix)[-1] <- paste0(data$metadata$SampleID, '_peak_area')
  
  if(save.table == FALSE){
    print('Table not saved')
  } else {
    write_csv(feature_abundance_matrix, save.table)
  }
  
  return(feature_abundance_matrix)
  
}

# -------------------------------------------------------------------------
extract_feature_definition <- function(data, feature_abundance_matrix, save.table = FALSE){
  ## get feature definitions and intensities
  feature_definition <- as.data.frame(featureDefinitions(data$spectra_centroid)) %>% 
    rownames_to_column(var = 'FeatureID') %>% 
    select(-peakidx) %>% 
    left_join(feature_abundance_matrix, by = 'FeatureID')
  
  if(save.table == FALSE){
    print('Table not saved')
  } else {
    write_csv(feature_definition, save.table)
  }
  
  return(feature_definition)
}

# -------------------------------------------------------------------------
extract_spectra_table <- function(data, save.table = FALSE){
  
  spectra_table <- fData(data$spectra_centroid) %>% 
    rownames_to_column(var = 'SpectraID')
  
  if(save.table == FALSE){
    print('Table not saved')
  } else {
    write_csv(spectra_table, save.table)
  }
  
  return(spectra_table)
}

# -------------------------------------------------------------------------
extract_MS2 <- function(data, save.spectra.consensus){
  # This function is based on a modified unction to write the feature name as the title of each MS2
  source('modified_writeMgfData.R')
  
  # This function uses custom functions for database compatibility from the Github Repo: 
  # https://github.com/jorainer/xcms-gnps-tools
  source("https://raw.githubusercontent.com/jorainer/xcms-gnps-tools/master/customFunctions.R")
  
  filteredMs2Spectra <- featureSpectra(data$spectra_centroid, return.type = "MSpectra", msLevel = 2)
  filteredMs2Spectra <- clean(filteredMs2Spectra, all = TRUE)
  filteredMs2Spectra <- formatSpectraForGNPS(filteredMs2Spectra) # this is one of the custom funtions
  filteredMs2Spectra_consensus <- combineSpectra(filteredMs2Spectra, 
                                                 fcol = "feature_id", 
                                                 method = consensusSpectrum, 
                                                 mzd = 0, 
                                                 minProp = 0.5, 
                                                 ppm = 25,
                                                 intensityFun = median,
                                                 mzFun = median)
  
  mod_writeMgfDataFile(filteredMs2Spectra_consensus, save.spectra.consensus)
  
  return(filteredMs2Spectra_consensus)
}


