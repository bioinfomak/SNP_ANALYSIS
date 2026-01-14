############################################################
## PUBLICATION-QUALITY 3-WAY VENN DIAGRAM (ggplot)
############################################################

library(ggVennDiagram)
library(ggplot2)

## -------------------------------
## 1. File paths
## -------------------------------
file1 <- "PATH"
file2 <- "PATH"
file3 <- "PATH"

## -------------------------------
## 2. Read & clean variant sets
## -------------------------------
v1 <- unique(na.omit(read.csv(file1, stringsAsFactors = FALSE)$variant_id))
v2 <- unique(na.omit(read.csv(file2, stringsAsFactors = FALSE)$variant_id))
v3 <- unique(na.omit(read.csv(file3, stringsAsFactors = FALSE)$variant_id))

stopifnot(length(v1) > 0, length(v2) > 0, length(v3) > 0)

## -------------------------------
## 3. Create Venn data
## -------------------------------
venn_list <- list(
  "M-SALS"   = v1,
  "F-SALS"   = v2,
  "Healthy"  = v3
)

## -------------------------------
## 4. Plot (publication-quality)
## -------------------------------
p <- ggVennDiagram(
  venn_list,
  label_alpha = 0,      # transparent labels background
  label_size = 6        # size of counts
) +
  scale_fill_gradient(
    low = "#F5F5F5",
    high = "#B22222"    # dark red for strong overlap
  ) +
  theme(
    plot.title = element_text(
      size = 18, face = "bold", hjust = 0.5
    ),
    legend.position = "none"
  ) +
  ggtitle("Variant Overlap Between M-SALS, F-SALS, and Healthy")

## -------------------------------
## 5. Save high-resolution figure
## -------------------------------
ggsave(
  filename = "Variant_Venn_ggplot_publication.png",
  plot = p,
  width = 7,
  height = 6,
  dpi = 300
)

## -------------------------------
## 6. Display plot
## -------------------------------
print(p)

# Read full CSVs
df_M_SALS <- read.csv(file1, stringsAsFactors = FALSE)
df_F_SALS <- read.csv(file2, stringsAsFactors = FALSE)
df_NORMAL <- read.csv(file3, stringsAsFactors = FALSE)

# Variant IDs
var1 <- unique(df_M_SALS$variant_id)
var2 <- unique(df_F_SALS$variant_id)
var3 <- unique(df_NORMAL$variant_id)

# Gene names
gene1 <- unique(df_M_SALS$gene_name)
gene2 <- unique(df_F_SALS$gene_name)
gene3 <- unique(df_NORMAL$gene_name)

common_variants_all <- Reduce(intersect, list(var1, var2, var3))

length(common_variants_all)
head(common_variants_all)

write.csv(
  data.frame(variant_id = common_variants_all),
  "Common_Variants_M_F_SALS_NORMAL.csv",
  row.names = FALSE
)

common_genes_all <- Reduce(intersect, list(gene1, gene2, gene3))

length(common_genes_all)
head(common_genes_all)

write.csv(
  data.frame(gene_name = common_genes_all),
  "Common_Genes_M_F_SALS_NORMAL.csv",
  row.names = FALSE
)

M_SALS_only <- setdiff(var1, union(var2, var3))
F_SALS_only <- setdiff(var2, union(var1, var3))
NORMAL_only <- setdiff(var3, union(var1, var2))

write.csv(data.frame(variant_id = M_SALS_only),
          "M_SALS_only_variants.csv", row.names = FALSE)

write.csv(data.frame(variant_id = F_SALS_only),
          "F_SALS_only_variants.csv", row.names = FALSE)

write.csv(data.frame(variant_id = NORMAL_only),
          "NORMAL_only_variants.csv", row.names = FALSE)

M_SALS_only_genes <- setdiff(gene1, union(gene2, gene3))
F_SALS_only_genes <- setdiff(gene2, union(gene1, gene3))
NORMAL_only_genes <- setdiff(gene3, union(gene1, gene2))

write.csv(data.frame(gene_name = M_SALS_only_genes),
          "M_SALS_only_genes.csv", row.names = FALSE)

write.csv(data.frame(gene_name = F_SALS_only_genes),
          "F_SALS_only_genes.csv", row.names = FALSE)

write.csv(data.frame(gene_name = NORMAL_only_genes),
          "NORMAL_only_genes.csv", row.names = FALSE)

