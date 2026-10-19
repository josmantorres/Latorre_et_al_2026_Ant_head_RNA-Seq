# ============================================================
# Differential expression analysis of Linepithema humile
# under carbohydrate restriction
#
# Comparison: Treated (low-sucrose diet) vs Control (normal diet)
# Statistical test: DESeq2 Wald test
# Log2 fold-change shrinkage: apeglm
# Low-expression filtering: HTSFilter
#
# Input files:
#   input/Raw_counts.txt
#   input/experimental_design.txt
#
# Output files:
#   results/DESeq2_Treated_vs_Control_apeglm.txt
#   results/Rlog_values.txt
# ============================================================


# 1. Load packages --------------------------------------------

suppressPackageStartupMessages({
  library(DESeq2)
  library(HTSFilter)
  library(apeglm)
})


# 2. Define input and output paths ----------------------------

input_dir <- "input"
output_dir <- "results"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

counts_file <- file.path(input_dir, "Raw_counts.txt")
metadata_file <- file.path(input_dir, "experimental_design.txt")


# 3. Import raw counts ----------------------------------------

counts_table <- read.delim(
  counts_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(counts_table) < 2) {
  stop("The count table must contain a gene ID column and sample columns.")
}

gene_ids <- counts_table[[1]]

if (anyNA(gene_ids) || any(gene_ids == "")) {
  stop("Missing or empty gene IDs were found in the count table.")
}

if (anyDuplicated(gene_ids)) {
  stop("Duplicated gene IDs were found in the count table.")
}

counts_matrix <- as.matrix(
  counts_table[, -1, drop = FALSE]
)

storage.mode(counts_matrix) <- "numeric"
rownames(counts_matrix) <- gene_ids

if (anyNA(counts_matrix) ||
    any(!is.finite(counts_matrix)) ||
    any(counts_matrix < 0) ||
    any(counts_matrix != round(counts_matrix))) {
  stop("Raw counts must be non-negative integer values without missing data.")
}

storage.mode(counts_matrix) <- "integer"


# 4. Import sample metadata -----------------------------------

metadata <- read.delim(
  metadata_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c("Sample", "Treatment")

if (!all(required_columns %in% colnames(metadata))) {
  stop("The metadata file must contain 'Sample' and 'Treatment' columns.")
}

if (anyNA(metadata$Sample) ||
    any(metadata$Sample == "") ||
    anyDuplicated(metadata$Sample)) {
  stop("Sample IDs must be present, non-empty, and unique.")
}

if (!setequal(colnames(counts_matrix), metadata$Sample)) {
  stop("Sample IDs do not match between the count matrix and metadata.")
}

# Reorder count columns to match the metadata.
counts_matrix <- counts_matrix[, metadata$Sample, drop = FALSE]

rownames(metadata) <- metadata$Sample

metadata$Treatment <- factor(
  metadata$Treatment,
  levels = c("Control", "Treated")
)

if (anyNA(metadata$Treatment)) {
  stop("Unexpected or missing treatment labels in the metadata.")
}

if (any(table(metadata$Treatment) == 0)) {
  stop("Both Control and Treated groups must contain samples.")
}

stopifnot(
  identical(colnames(counts_matrix), rownames(metadata))
)


# 5. Create the DESeq2 dataset -------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = metadata,
  design = ~ Treatment
)

# Set Control as the reference group.
dds$Treatment <- relevel(dds$Treatment, ref = "Control")


# 6. Generate rlog-transformed expression values --------------

# The transformation is applied before HTSFilter, as in the
# original analysis. The output therefore includes all input genes.
rlog_data <- rlog(
  dds,
  blind = TRUE
)

write.table(
  assay(rlog_data),
  file = file.path(output_dir, "Rlog_values.txt"),
  sep = "\t",
  dec = ".",
  quote = FALSE,
  col.names = NA
)


# 7. Filter low-expression genes using HTSFilter ---------------

hts_result <- HTSFilter(
  dds,
  s.len = 100,
  plot = FALSE
)

dds_filtered <- hts_result$filteredData

if (nrow(dds_filtered) == 0) {
  stop("HTSFilter removed all genes. Check the filtering settings.")
}


# 8. Differential expression analysis ------------------------

dds_filtered <- DESeq(
  dds_filtered,
  test = "Wald"
)

# Verify the expected coefficient.
coef_name <- "Treatment_Treated_vs_Control"

if (!coef_name %in% resultsNames(dds_filtered)) {
  stop(
    "The expected coefficient was not found. Available coefficients: ",
    paste(resultsNames(dds_filtered), collapse = ", ")
  )
}


# 9. Extract results and shrink log2 fold changes --------------

results_raw <- results(
  dds_filtered,
  name = coef_name,
  independentFiltering = FALSE,
  cooksCutoff = TRUE
)

results_apeglm <- lfcShrink(
  dds_filtered,
  coef = coef_name,
  res = results_raw,
  type = "apeglm"
)

results_table <- as.data.frame(results_apeglm)
results_table$GeneID <- rownames(results_table)

results_table <- results_table[
  ,
  c("GeneID", setdiff(colnames(results_table), "GeneID")),
  drop = FALSE
]

write.table(
  results_table,
  file = file.path(
    output_dir,
    "DESeq2_Treated_vs_Control_apeglm.txt"
  ),
  sep = "\t",
  dec = ".",
  row.names = FALSE,
  quote = FALSE
)


# 10. Report completion ---------------------------------------

cat("Differential expression analysis completed successfully.\n")
cat("Output files:\n")
cat(" - results/DESeq2_Treated_vs_Control_apeglm.txt\n")
cat(" - results/Rlog_values.txt\n")
