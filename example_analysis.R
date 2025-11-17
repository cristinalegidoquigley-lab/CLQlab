# Example Analysis Using Synthetic MicrobLiver Data
# This script demonstrates how to use the synthetic datasets with your analysis functions

# ==============================================================================
# SETUP
# ==============================================================================

# Set your working directory to where you have the files
# setwd("path/to/your/folder")

# Load the analysis functions
source("functions.R")

# Define directories
data_dir <- "preprocessData/"
results_dir <- "results/"
models_dir <- "models/"

# Create output directories if they don't exist
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(models_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# EXAMPLE 1: Load and Explore a Dataset
# ==============================================================================

cat("\n=== EXAMPLE 1: Loading Data ===\n")

# Load metabolites data (hepatic samples only)
metabolites_hepatic <- readFile(
  name = "plasma_metabolites",
  data_dir = data_dir,
  sampleSite = "hepatic"
)

cat("Loaded hepatic metabolites:\n")
cat("  - Samples:", nrow(metabolites_hepatic), "\n")
cat("  - Features:", ncol(metabolites_hepatic) - 7, "\n")
cat("  - Phenotype breakdown:\n")
print(table(metabolites_hepatic$pheno))

# Load systemic samples
metabolites_systemic <- readFile(
  name = "plasma_metabolites",
  data_dir = data_dir,
  sampleSite = "systemic"
)

cat("\nLoaded systemic metabolites:\n")
cat("  - Samples:", nrow(metabolites_systemic), "\n")

# ==============================================================================
# EXAMPLE 2: Run Linear Mixed Models
# ==============================================================================

cat("\n=== EXAMPLE 2: Linear Mixed Models ===\n")

# Analyze metabolites over time in hepatic samples
# This tests for metabolite changes at 60 and 180 min vs baseline (0 min)
# Separate analysis for each phenotype group (Healthy, ALD, NAFLD)

cat("\nRunning LMM on hepatic metabolites...\n")
results_metabolites_hepatic <- LMM_time(
  name = "plasma_metabolites",
  ID = "PatientID",
  fixedEffect = "time",
  groupBy = "pheno",
  data_dir = data_dir,
  results_dir = results_dir,
  models_dir = models_dir,
  sampleSite = "hepatic",
  overwrite = TRUE,
  padj_method = "fdr"
)

cat("Analysis complete!\n")
cat("  - Total results:", nrow(results_metabolites_hepatic), "\n")
cat("  - Significant (padj<0.05):", sum(results_metabolites_hepatic$padj < 0.05, na.rm=TRUE), "\n")

# Show the most significant results
sig_results <- results_metabolites_hepatic[
  results_metabolites_hepatic$padj < 0.05 & 
  !grepl("Intr", results_metabolites_hepatic$fixed_effect), 
]

if(nrow(sig_results) > 0) {
  cat("\nTop 5 significant changes:\n")
  sig_results <- sig_results[order(sig_results$padj), ]
  print(head(sig_results[, c("molecule", "fixed_effect", "Estimate", "pvalue", "padj")], 5))
}

# ==============================================================================
# EXAMPLE 3: Analyze Multiple Datasets
# ==============================================================================

cat("\n=== EXAMPLE 3: Analyzing All Omics Layers ===\n")

# You can run the same analysis on different omics layers
omics_types <- c("plasma_metabolites", "plasma_lipidomics", "plasma_proteomics")
sample_site <- "hepatic"

for(omics_type in omics_types) {
  cat("\nAnalyzing", omics_type, "...\n")
  
  results <- LMM_time(
    name = omics_type,
    ID = "PatientID",
    fixedEffect = "time",
    groupBy = "pheno",
    data_dir = data_dir,
    results_dir = results_dir,
    models_dir = models_dir,
    sampleSite = sample_site,
    overwrite = TRUE,
    padj_method = "fdr"
  )
  
  n_sig <- sum(results$padj < 0.05, na.rm=TRUE)
  cat("  - Significant features:", n_sig, "\n")
}

# ==============================================================================
# EXAMPLE 4: Create Visualizations
# ==============================================================================

cat("\n=== EXAMPLE 4: Creating Volcano Plots ===\n")

# Create volcano plots for the metabolites results
# This requires ggplot2 and ggrepel packages

# Uncomment to run (requires packages):
# plots <- volcanoPlot_time(results_metabolites_hepatic, pvalue = "padj")
# 
# # Save plots
# for(plot_name in names(plots)) {
#   ggsave(
#     filename = paste0("volcano_", plot_name, ".png"),
#     plot = plots[[plot_name]],
#     width = 8,
#     height = 6,
#     dpi = 300
#   )
# }

cat("\nVolcano plot function available - requires ggplot2 and ggrepel packages\n")

# ==============================================================================
# EXAMPLE 5: Compare Hepatic vs Systemic
# ==============================================================================

cat("\n=== EXAMPLE 5: Comparing Sampling Sites ===\n")

# Run analysis on systemic samples
results_metabolites_systemic <- LMM_time(
  name = "plasma_metabolites",
  ID = "PatientID",
  fixedEffect = "time",
  groupBy = "pheno",
  data_dir = data_dir,
  results_dir = results_dir,
  models_dir = models_dir,
  sampleSite = "systemic",
  overwrite = TRUE,
  padj_method = "fdr"
)

# Compare number of significant features
n_sig_hepatic <- sum(results_metabolites_hepatic$padj < 0.05, na.rm=TRUE)
n_sig_systemic <- sum(results_metabolites_systemic$padj < 0.05, na.rm=TRUE)

cat("\nComparison of sampling sites:\n")
cat("  - Hepatic significant features:", n_sig_hepatic, "\n")
cat("  - Systemic significant features:", n_sig_systemic, "\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("\nGenerated files:\n")
cat("  - Results:", list.files(results_dir, pattern=".rds"), "\n")
cat("  - Models:", list.files(models_dir, pattern=".rds"), "\n")

cat("\nNext steps:\n")
cat("  1. Explore the results objects\n")
cat("  2. Create volcano plots and other visualizations\n")
cat("  3. Integrate multiple omics layers\n")
cat("  4. Create chord diagrams for multi-omics correlations\n")
cat("  5. Generate circular heatmaps\n")

cat("\nFor more details, see the README.md in preprocessData/\n")
