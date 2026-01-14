############################################################
## COMMON VARIANT FINDER WITH GENE SYMBOL
## Input : Folder with annotated.csv files
## Output: common_variant_annotated.csv
############################################################

## -------------------------------
## 1. Set input directory
## -------------------------------
input_dir <- "PATH"

## -------------------------------
## 2. List annotated.csv files
## -------------------------------
csv_files <- list.files(
  path = input_dir,
  pattern = "annotated\\.csv$",
  full.names = TRUE
)

if (length(csv_files) < 2) {
  stop("At least two annotated.csv files are required.")
}

cat("Total annotated.csv files found:", length(csv_files), "\n")

## -------------------------------
## 3. Read all annotated files
## -------------------------------
annotated_list <- lapply(csv_files, function(file) {
  df <- read.csv(file, stringsAsFactors = FALSE)
  
  required_cols <- c("variant_id", "gene_name")   # <- change if needed
  
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0) {
    stop(paste("Missing columns in", basename(file), ":", paste(missing, collapse = ", ")))
  }
  
  df
})

## -------------------------------
## 4. Find common variant_id
## -------------------------------
variant_lists <- lapply(annotated_list, function(df) unique(df$variant_id))
common_variant <- Reduce(intersect, variant_lists)

cat("Common variants found:", length(common_variant), "\n")

## -------------------------------
## 5. Keep annotation (gene symbol)
##    Use first file as reference
## -------------------------------
final_annotated <- annotated_list[[1]][
  annotated_list[[1]]$variant_id %in% common_variant,
  c("variant_id", "gene_name")
]

## Remove duplicates (safety)
final_annotated <- unique(final_annotated)

## -------------------------------
## 6. Save output
## -------------------------------
write.csv(
  final_annotated,
  file = "common_variant_annotated.csv",
  row.names = FALSE,
  quote = FALSE
)

## -------------------------------
## 7. Summary
## -------------------------------
cat("=====================================\n")
cat("Final variants:", nrow(final_annotated), "\n")
cat("Output: common_variant_annotated.csv\n")
cat("=====================================\n")

head(final_annotated, 20)

############################################################
## UPSET PLOT FOR COMMON VARIANTS
## Input : Multiple annotated.csv files
############################################################

## -------------------------------
## 1. Load libraries
## -------------------------------
suppressPackageStartupMessages({
  library(UpSetR)
  library(tools)
})

## -------------------------------
## 2. Set input directory
## -------------------------------
input_dir <- "PATH"

## -------------------------------
## 3. List annotated.csv files
## -------------------------------
csv_files <- list.files(
  input_dir,
  pattern = "annotated\\.csv$",
  full.names = TRUE
)

stopifnot(length(csv_files) >= 2)

sample_names <- file_path_sans_ext(basename(csv_files))

## -------------------------------
## 4. Read variant_id lists
## -------------------------------
variant_list <- lapply(csv_files, function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  unique(df$variant_id)
})

names(variant_list) <- sample_names

## -------------------------------
## 5. Create UpSet input
## -------------------------------
upset_input <- fromList(variant_list)

## -------------------------------
## 6. Plot UpSet
## -------------------------------
upset(
  upset_input,
  sets = sample_names,
  order.by = "freq",
  main.bar.color = "black",
  sets.bar.color = "gray40",
  text.scale = 1.5
)
