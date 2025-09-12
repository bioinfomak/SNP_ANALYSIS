# Load required libraries
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(SNPlocs.Hsapiens.dbSNP156.GRCh38)
library(org.Hs.eg.db)
library(GenomicRanges)
library(readr)
library(dplyr)

# ==== Set your VCF folder path here ====
vcf_folder <- "PATH/FILE"  # <<< CHANGE THIS to your folder path

# List all .vcf files in the folder (case-insensitive)
vcf_files <- list.files(path = vcf_folder, pattern = "\\.vcf$", full.names = TRUE, ignore.case = TRUE)

# Create output folder for CSV results
output_folder <- file.path(vcf_folder, "GeneOverlapResults")
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# Load gene annotation once
gr_genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Map ENTREZ IDs to gene symbols once
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = as.character(gr_genes$gene_id),
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")
gr_genes$gene_name <- gene_symbols

# Process each VCF file
for (vcf_file in vcf_files) {
  cat("Processing:", vcf_file, "\n")
  
  # Read VCF
  vcf <- readVcf(vcf_file, genome = "hg38")
  gr_variants <- rowRanges(vcf)
  
  # Find overlaps with genes
  hits <- findOverlaps(gr_variants, gr_genes, ignore.strand = TRUE)
  
  # Collapse ALT alleles for variants with overlaps
  alt_alleles <- sapply(mcols(gr_variants)$ALT[queryHits(hits)], function(x) {
    paste(as.character(x), collapse = ",")
  })
  
  # Create overlap data frame
  overlap_table <- data.frame(
    variant_id = names(gr_variants)[queryHits(hits)],
    chr        = as.character(seqnames(gr_variants)[queryHits(hits)]),
    pos        = start(gr_variants)[queryHits(hits)],
    ref        = as.character(mcols(gr_variants)$REF[queryHits(hits)]),
    alt        = alt_alleles,
    gene_id    = gr_genes$gene_id[subjectHits(hits)],
    gene_name  = gr_genes$gene_name[subjectHits(hits)],
    stringsAsFactors = FALSE
  )
  
  # Build output file name and path
  output_file <- sub("\\.vcf$", "_VCF_Gene.csv", basename(vcf_file), ignore.case = TRUE)
  output_path <- file.path(output_folder, output_file)
  
  # Write CSV file
  write.csv(overlap_table, file = output_path, row.names = FALSE)
  
  cat("Saved results to:", output_path, "\n\n")
}

cat("All files processed!\n")


# FIND SNP AND INDEL FROM VCF_GENE_FILE

# Step 1: Read your CSV file
df <- read_csv("VCF_FILE")

# Step 2: Check column names
print(colnames(df))
stopifnot(all(c("chr", "pos", "ref", "alt") %in% colnames(df)))

# Step 3: Normalize chromosome names
df$chr <- sub("^chr", "", df$chr)

# Step 4: Flag indels (insertions or deletions)
df <- df %>%
  mutate(
    ref_length = nchar(ref),
    alt_length = nchar(alt),
    is_indel = (ref_length != 1 | alt_length != 1)
  )

# Step 5: Create GRanges object for SNP positions only (non-indels)
gr <- GRanges(
  seqnames = df$chr[!df$is_indel],
  ranges = IRanges(start = df$pos[!df$is_indel], end = df$pos[!df$is_indel])
)

# Step 6: Load dbSNP SNP-only database
snp_db <- SNPlocs.Hsapiens.dbSNP156.GRCh38

# Step 7: Sync seqlevels to avoid warnings
common_seqlevels <- intersect(seqlevels(gr), seqlevels(snp_db))
if (length(common_seqlevels) == 0) {
  stop("No common chromosomes found between your data and dbSNP database!")
}
gr <- keepSeqlevels(gr, common_seqlevels, pruning.mode = "coarse")

# Step 8: Query SNPs overlapping your positions
snp_hits <- snpsByOverlaps(snp_db, gr)
snp_df <- as.data.frame(snp_hits)

# Step 9: Aggregate SNP rsIDs per position (optional if multiple rsIDs exist)
snp_agg <- snp_df %>%
  group_by(seqnames, pos) %>%
  summarise(rsIDs = paste(unique(RefSNP_id), collapse = ";")) %>%
  ungroup()

# Step 10: Join SNP rsIDs back to original dataframe
# SNP rsIDs will only be matched for non-indels
result <- df %>%
  left_join(snp_agg, by = c("chr" = "seqnames", "pos" = "pos"))

# Step 11: Create a note column for indels
result <- result %>%
  mutate(
    rsIDs = ifelse(is_indel, NA, rsIDs),  # remove SNP rsIDs from indels just in case
    indel_note = ifelse(is_indel,
                        "Indel not annotated (SNP database only)",
                        "SNP rsID matched if available")
  )

# Step 12: Save output to CSV
write.csv(result, "SS4_155_rsid_matched_with_indel_flag.csv", row.names = FALSE)

# Optional: View summary
print(head(result))
table(result$is_indel)
