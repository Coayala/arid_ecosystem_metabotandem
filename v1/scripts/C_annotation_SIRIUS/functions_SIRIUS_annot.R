################################
# FUNCTIONS SIRIUS ANNOTATION #
###############################

# -------------------------------------------------------------------------
separate_MS2 <- function(mgf_file, outdir, n){
  mgf <- read_lines(mgf_file)
  
  if(!dir.exists(outdir)){
    dir.create(outdir)
  }
  
  accum <- c()
  
  start <- 'FALSE'
  stop <-  'FALSE'
  group <- 1
  counter <- 0
  
  for(i in 1:length(mgf)){
    if(str_detect(mgf[i], 'BEGIN')){
      start <- "TRUE"
      stop <-  'FALSE'
      counter <- counter + 1
    } 
    
    if(str_detect(mgf[i], 'END')){
      start <- "FALSE"
      stop <- "TRUE"
    }
    
    if(start == 'TRUE'){
      accum <- c(accum, mgf[i])
    } 
    
    if(stop == 'TRUE'){
      accum <- c(accum, mgf[i])
      if(counter == n){
        counter <- 0
        write_lines(accum, paste0(outdir, 'group_', group, '.mgf'))
        group <- group + 1
        accum <- c()
      }
    }
    if(i == length(mgf)){
      write_lines(accum, paste0(outdir, 'group_', group, '.mgf'))
    }
  }
  
}

# -------------------------------------------------------------------------
annotate_SIRIUS <- function(groupdir, outdir){
  
  mgf_list <- list.files(groupdir, pattern = '*.mgf', full.names = TRUE)
  
  
  
  for(file in mgf_list){
    
    file <- mgf_list[1]
    name <- unlist(str_split(file, '\\/'))
    name <- str_remove(name[length(name)], '.mgf')
    outfile <- file.path(outdir, name)
    system2('./run_sirius_cli.sh', outfile, file)
  }
}

# -------------------------------------------------------------------------
merge_SIRIUS_files <- function(annot_dir){
  
  id_files <- list.files(annot_dir, pattern = 'compound_identifications.tsv', recursive = TRUE, full.names = TRUE)
  canopus_files <- list.files(annot_dir, pattern = 'canopus_summary.tsv', recursive = TRUE, full.names = TRUE)
  
  id_merged <- as.data.frame(matrix(ncol = 6, nrow = 0))
  colnames(id_merged) <- c('FeatureID', 'molecularFormula', 'adduct', 'InChI', 'smiles', 'links')
  
  for(file in id_files){
    temp_df <- read_tsv(file) %>% 
      mutate(FeatureID = str_extract(id, 'FT.*')) %>% 
      select(FeatureID, molecularFormula, adduct, InChI, smiles, links)
    
    id_merged <- rbind(id_merged, temp_df)
  }
  
  canopus_merged <- as.data.frame(matrix(ncol = 5, nrow = 0))
  colnames(canopus_merged) <- c('FeatureID', "all classifications", 'superclass', 'class', 'subclass')
  
  
  for(file in canopus_files){
    temp_df <- read_tsv(file) %>% 
      mutate(FeatureID = str_extract(name, 'FT.*')) %>% 
      select(FeatureID, `all classifications`, superclass, class, subclass)
    
    canopus_merged <- rbind(canopus_merged, temp_df)
  }
  
  annotation_df <- left_join(id_merged, canopus_merged, by = 'FeatureID')
  
  return(annotation_df)
}

# -------------------------------------------------------------------------
fix_entries <- function(annotation_df, save.table = FALSE){
  fixed_df <- annotation_df %>% 
    mutate(HMDB = str_extract(links, 'HMDB:\\([0-9| ]+\\)'),
           HMDB = str_trim(str_remove_all(HMDB, 'HMDB:\\(|\\)'), side = 'both'),
           YMDB = str_extract(links, 'YMDB:\\([0-9| ]+\\)'),
           YMDB = str_remove_all(YMDB, 'YMDB:\\(|\\)'),
           KNApSAcK = str_extract(links, 'KNApSAcK:\\([0-9| ]+\\)'),
           KNApSAcK = str_remove_all(KNApSAcK, 'KNApSAcK:\\(|\\)'),
           CHEBI = str_extract(links, 'CHEBI:\\([0-9| ]+\\)'),
           CHEBI = str_remove_all(CHEBI, 'CHEBI:\\(|\\)'),
           PlantCyc = str_extract(links, 'Plantcyc:\\([A-Z|0-9| |-]+\\)'),
           PlantCyc = str_remove_all(PlantCyc, 'Plantcyc:\\(|\\)'),
           PlantCyc = str_replace_all(PlantCyc, ' ', ','),
           BioCyc = str_extract(links, 'Biocyc:\\([A-Z|0-9| |-]+\\)'),
           BioCyc = str_remove_all(BioCyc, 'Biocyc:\\(|\\)'),
           BioCyc = str_replace_all(BioCyc, ' ', ','),
           KEGG = str_extract(links, 'KEGG:\\(C[0-9| ]+\\)'),
           KEGG = str_remove_all(KEGG, 'KEGG:\\(|\\)'),
           KEGG = str_replace_all(KEGG, ' ', ','),
           COCONUT = str_extract(links, 'COCONUT:\\(CNP[0-9| ]+\\)'),
           COCONUT = str_remove_all(COCONUT, 'COCONUT:\\(|\\)'),
           COCONUT = str_replace_all(COCONUT, ' ', ','),
           PubChem_CID = str_extract(links, 'PubChem:\\([0-9| ]+\\)'),
           PubChem_CID = str_remove_all(PubChem_CID, 'PubChem:\\(|\\)'),
           PubChem_CID = str_replace_all(PubChem_CID, ' ', ','))
  
  for(i in 1:nrow(fixed_df)){
    if(!is.na(fixed_df$HMDB[i])){
      fixed_df$HMDB[i] <- paste(sprintf('HMDB%07d', as.numeric(unlist(str_split(fixed_df$HMDB[i], ' ')))), collapse = ',')
    }
    if(!is.na(fixed_df$YMDB[i])){
      fixed_df$YMDB[i] <- paste(sprintf('YMDB%05d', as.numeric(unlist(str_split(fixed_df$YMDB[i], ' ')))), collapse = ',')
    }
    if(!is.na(fixed_df$KNApSAcK[i])){
      fixed_df$KNApSAcK[i] <- paste(sprintf('C%08d', as.numeric(unlist(str_split(fixed_df$KNApSAcK[i], ' ')))), collapse = ',')
    }
    if(!is.na(fixed_df$CHEBI[i])){
      fixed_df$CHEBI[i] <- paste(paste0('CHEBI:', as.numeric(unlist(str_split(fixed_df$CHEBI[i], ' ')))), collapse = ',')
    }
  }
  
  if(save.table == FALSE){
    print('Table not saved')
  } else {
    write_csv(fixed_df, save.table)
  }
  
  return(fixed_df)
} 

# -------------------------------------------------------------------------
explode_table <- function(annotation_df, explode_by = 'PubChem_CID'){
  
 exploded_df <- annotation_df %>% 
   separate_rows('PubChem_CID', sep = ',')
 
 return(exploded_df)
}

# -------------------------------------------------------------------------
retrieve_pubchem_hits <- function(annotation_df, save.table = FALSE){
  cid_hits <- unique(annotation_df$PubChem_CID)
  cid_hits <- cid_hits[!is.na(cid_hits)]
  
  
  cid_df <- tibble(CID = NA, MolecularFormula = NA, ExactMass = NA, InChKey = NA, .rows = 0)
  
  for(i in 1:length(cid_hits)){
   
    re <- Get(paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/', i, '/property/MolecularFormula/CSV'))
    mf = str_remove_all(unlist(str_split(re$content, ','))[3], '\\\n|\\"|\\\\')
    
    re <- Get(paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/', i, '/property/ExactMass/CSV'))
    em = str_remove_all(unlist(str_split(re$content, ','))[3], '\\\n|\\"|\\\\')
    
    re <- Get(paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/', i, '/property/InChIKey/CSV'))
    ik = str_remove_all(unlist(str_split(re$content, ','))[3], '\\\n|\\"|\\\\')
    
    
    temp_df <- tibble(CID = cid_hits[i], 
                      MolecularFormula = mf, 
                      ExactMass = em, 
                      InChKey = ik)
    
    cid_df <- rbind(cid_df, temp_df)
    
    print(paste0('[', i, '/', length(cid_hits), '] Retrieving information for CID ', cid_hits[i]))
    
    if(i %% 10 == 0){
      Sys.sleep(5)
    }
  }
  
  final_df <- left_join(annotation_df, cid_df, by = c('PubChem_CID' = 'CID'))
  
  if(save.table == FALSE){
    print('Table not saved')
  } else {
    write_csv(final_df, save.table)
  }
  
  return(cid_df)
  
}


