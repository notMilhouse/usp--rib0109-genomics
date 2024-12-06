# Processing RNAseq Data

# Define directories for jovem and naojovem RNAseq data
rnaseq_dirs <- list.dirs("data/rnaseq", recursive = FALSE)
rnaseq_categories <- c("jovem", "naojovem")

# Initialize a data frame for RNAseq data
rnaseq_data <- data.frame(gene_id = character())

# Loop through jovem and naojovem categories for RNAseq
for (category in rnaseq_categories) {
  cat_dir <- rnaseq_dirs[grepl(category, rnaseq_dirs)]
  
  # Ensure category directory exists
  if (length(cat_dir) == 0) {
    cat(paste("No directory found for category:", category, "\n"))
    next
  }
  
  # Get list of .tsv files for the category
  files <- list.files(cat_dir, pattern = "\\.tsv$", full.names = TRUE)
  
  # Sum counts across all samples for this category
  category_counts <- NULL
  for (file in files) {
    tmp_data <- read.delim(file, skip = 1)[-c(1:4), c(1, 4)]  # Read gene_id and counts
    colnames(tmp_data) <- c("gene_id", "count")
    tmp_data$count <- as.numeric(tmp_data$count)
    
    if (is.null(category_counts)) {
      category_counts <- tmp_data
    } else {
      category_counts <- merge(category_counts, tmp_data, by = "gene_id", all = TRUE)
      category_counts$count <- rowSums(category_counts[, -1], na.rm = TRUE)
      category_counts <- category_counts[, c("gene_id", "count")]
    }
  }
  
  # Rename count column to reflect the category
  colnames(category_counts)[2] <- paste0(category, "_count")
  
  # Merge category counts with the main RNAseq data
  if (ncol(rnaseq_data) == 1) {  # Only gene_id column present
    rnaseq_data <- category_counts
  } else {
    rnaseq_data <- merge(rnaseq_data, category_counts, by = "gene_id", all = TRUE)
  }
}

# Save the final RNAseq data
saveRDS(rnaseq_data, file = "data/processed_data/rnaseq_data.rds")
print(head(rnaseq_data))

# Processing MAF Data

# Define directories for jovem and naojovem MAF data
maf_dirs <- list.dirs("data/maf", recursive = FALSE)
maf_categories <- c("jovem", "naojovem")

# Initialize a data frame for MAF data
maf_data <- NULL

# Loop through jovem and naojovem categories for MAF
for (category in maf_categories) {
  cat_dir <- maf_dirs[grepl(category, maf_dirs)]
  
  # Ensure category directory exists
  if (length(cat_dir) == 0) {
    cat(paste("No directory found for category:", category, "\n"))
    next
  }
  
  # Get list of .maf.gz files for the category
  files <- list.files(cat_dir, pattern = "maf\\.gz$", full.names = TRUE)
  
  # Initialize category-specific data frame
  category_data <- NULL
  
  for (file in files) {
    # Read data, skip first 7 rows and select relevant columns
    tmp_data <- read.delim(file, skip = 7)[, c(1, 9, 16)]
    colnames(tmp_data) <- c("Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode")
    
    # Add category as a new column
    tmp_data$Category <- category
    
    if (is.null(category_data)) {
      category_data <- tmp_data
    } else {
      category_data <- rbind(category_data, tmp_data)
    }
  }
  
  # Merge category data into the main maf_data
  if (is.null(maf_data)) {
    maf_data <- category_data
  } else {
    maf_data <- rbind(maf_data, category_data)
  }
}

# Filter MAF data for relevant mutations
muts <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Missense_Mutation", "Nonsense_Mutation")
maf_data_filtered <- maf_data[maf_data$Variant_Classification %in% muts, ]

# Save the final MAF data
saveRDS(maf_data, file = "data/processed_data/maf_data.rds")
saveRDS(maf_data_filtered, file = "data/processed_data/maf_data_filtered.rds")

# Final Output Preview
cat("RNAseq Data Preview:\n")
print(head(rnaseq_data))
cat("\nMAF Data Preview:\n")
print(head(maf_data))
