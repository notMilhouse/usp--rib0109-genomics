# Define helper function for RNAseq data processing
load_rnaseq_data <- function(base_path, category) {
  category_path <- file.path(base_path, "rnaseq", category)
  rnaseq_dirs <- list.files(category_path)
  
  if (length(rnaseq_dirs) == 0) {
    stop(paste("No directories found in", category_path))
  }
  
  # Read and combine RNAseq data
  raw_counts <- lapply(rnaseq_dirs, function(dir) {
    file_path <- list.files(file.path(category_path, dir), pattern = '\\.tsv$', full.names = TRUE)
    if (length(file_path) == 0) {
      stop(paste("No TSV files found in", file.path(category_path, dir)))
    }
    tmp <- read.delim(file_path, skip = 1)[-c(1:4), c(1, 4)]
    colnames(tmp) <- c("gene_id", dir)
    tmp
  })
  
  # Merge data frames by gene_id
  raw_counts <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), raw_counts)
  
  # Convert to matrix and set rownames
  rownames(raw_counts) <- raw_counts$gene_id
  as.matrix(raw_counts[, -1])
}

# Define helper function for MAF data processing
load_maf_data <- function(base_path, category) {
  category_path <- file.path(base_path, "maf", category)
  maf_dirs <- list.files(category_path)
  
  if (length(maf_dirs) == 0) {
    stop(paste("No directories found in", category_path))
  }
  
  # Read and combine MAF data
  maf <- do.call(rbind, lapply(maf_dirs, function(dir) {
    file_path <- list.files(file.path(category_path, dir), pattern = 'maf.gz$', full.names = TRUE)
    if (length(file_path) == 0) {
      stop(paste("No MAF files found in", file.path(category_path, dir)))
    }
    tmp <- read.delim(file_path, skip = 7)[, c(1, 9, 16)]
    colnames(tmp) <- c("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification")
    tmp
  }))
  
  maf
}

# Main function to process data
process_data <- function(base_path, categories, data_type) {
  processed_data <- lapply(categories, function(category) {
    if (data_type == "rnaseq") {
      load_rnaseq_data(base_path, category)
    } else if (data_type == "maf") {
      load_maf_data(base_path, category)
    }
  })
  
  # Save processed data for each category
  names(processed_data) <- categories
  for (category in categories) {
    saveRDS(processed_data[[category]], file = file.path(base_path, "processed_data", paste0("TCGA_UCS_", toupper(data_type), "_", toupper(category), ".rds")))
  }
  
  processed_data
}

# Paths and categories
base_path <- "data"
categories <- c("jovem", "naojovem")

# Process RNAseq data
rnaseq_data <- process_data(base_path, categories, "rnaseq")

# Process MAF data
maf_data <- process_data(base_path, categories, "maf")

# Optional: Combine RNAseq data for integration
combined_rnaseq <- do.call(cbind, rnaseq_data)
saveRDS(combined_rnaseq, file = file.path(base_path, "processed_data", "TCGA_UCS_RNASEQ_COMBINED.rds"))
