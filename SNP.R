library(readr)
library(GenomicRanges)
library(SNPlocs.Hsapiens.dbSNP156.GRCh38)
library(dplyr)
library(tools)  # for file_path_sans_ext()

# ---- USER SETTINGS ----
# Set the folder containing your .csv files
input_folder <- "PATH"
output_folder <- input_folder  # You can change this if you want outputs elsewhere

# List all CSV files in the folder
csv_files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)

# Load dbSNP database (only once)
snp_db <- SNPlocs.Hsapiens.dbSNP156.GRCh38

# Function to process each CSV file
process_snp_file <- function(file_path) {
  cat("Processing:", file_path, "\n")
  
  # Step 1: Read input CSV
  df <- read.csv(file_path)
  
  # Step 2: Check required columns
  required_cols <- c("chr", "pos", "ref", "alt")
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0) {
    warning(paste("Skipping file due to missing columns:", paste(missing, collapse = ", ")))
    return(NULL)
  }
  
  # Step 3: Normalize chromosome names (remove "chr" if present)
  df$chr <- sub("^chr", "", df$chr)
  
  # Step 4: Flag indels
  df <- df %>%
    mutate(
      ref_length = nchar(ref),
      alt_length = nchar(alt),
      is_indel = (ref_length != 1 | alt_length != 1)
    )
  
  # Step 5: Create GRanges for SNPs only
  snp_only_df <- df[!df$is_indel, ]
  if (nrow(snp_only_df) == 0) {
    cat("No SNPs found in this file. Skipping rsID annotation.\n")
    df$rsIDs <- NA
    df$indel_note <- "Indel not annotated (SNP database only)"
  } else {
    gr <- GRanges(
      seqnames = snp_only_df$chr,
      ranges = IRanges(start = snp_only_df$pos, end = snp_only_df$pos)
    )
    
    # Step 6: Match seqlevels
    common_seqlevels <- intersect(seqlevels(gr), seqlevels(snp_db))
    if (length(common_seqlevels) == 0) {
      warning("No common chromosomes found with dbSNP. Skipping.")
      return(NULL)
    }
    gr <- keepSeqlevels(gr, common_seqlevels, pruning.mode = "coarse")
    
    # Step 7: Query SNPs by overlap
    gr <- sort(gr)  # optional but good practice
    snp_hits <- snpsByOverlaps(snp_db, gr)
    snp_df <- as.data.frame(snp_hits)
    
    # Step 8: Aggregate rsIDs by position
    snp_agg <- snp_df %>%
      group_by(seqnames, pos) %>%
      summarise(rsIDs = paste(unique(RefSNP_id), collapse = ";"), .groups = "drop")
    
    # Step 9: Join rsIDs (handle case where snp_agg might be empty)
    if (nrow(snp_agg) == 0) {
      df$rsIDs <- NA
      warning(paste("No SNP rsIDs matched in file:", basename(file_path)))
    } else {
      df <- df %>%
        left_join(snp_agg, by = c("chr" = "seqnames", "pos" = "pos"))
    }
    
    # Step 10: Add annotation notes
    df <- df %>%
      mutate(
        rsIDs = ifelse(is_indel, NA, rsIDs),  # Remove rsIDs from indels
        indel_note = ifelse(is_indel,
                            "Indel not annotated (SNP database only)",
                            "SNP rsID matched if available")
      )
  }
  
  # Step 11: Save output file
  output_file <- file.path(output_folder,
                           paste0(file_path_sans_ext(basename(file_path)), "_annotated.csv"))
  write.csv(df, output_file, row.names = FALSE)
  cat("Saved to:", output_file, "\n\n")
}

# ---- BATCH PROCESSING ----
cat("Starting batch SNP annotation...\n")
lapply(csv_files, process_snp_file)
cat("Batch processing complete.\n")
