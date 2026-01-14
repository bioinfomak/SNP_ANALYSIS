############################################################
## VENN DIAGRAM FOR 3 CSV FILES (VARIANT_ID) – WITH COUNTS
############################################################

library(VennDiagram)
library(grid)

## -------------------------------
## 1. File paths
## -------------------------------
file1 <- "PATH"
file2 <- "PATH"
file3 <- "PATH"

## -------------------------------
## 2. Read variant IDs
## -------------------------------
var1 <- unique(read.csv(file1, stringsAsFactors = FALSE)$variant_id)
var2 <- unique(read.csv(file2, stringsAsFactors = FALSE)$variant_id)
var3 <- unique(read.csv(file3, stringsAsFactors = FALSE)$variant_id)

## -------------------------------
## 3. Create Venn plot object
## -------------------------------
venn.plot <- venn.diagram(
  x = list(
    M_SALS = var1,
    F_SALS = var2,
    NORMAL = var3
  ),
  filename = NULL,
  fill = c("orange", "purple", "cyan"),
  alpha = 1,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = c(-20, 20, -180),  # adjusted third label position
  cat.dist = c(0.05, 0.05, 0.05),  # adjusted third label distance
  main = "Variant Overlap Between Groups"
)

png("Variant_Venn_3samples_with_counts.png",
    width = 1500, height = 1500, res = 1200)
grid.draw(venn.plot)
dev.off()

grid.draw(venn.plot)

