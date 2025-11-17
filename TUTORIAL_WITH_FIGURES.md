## 📚 Table of Contents
1. [Introduction](#introduction)
2. [Study Design](#study-design)
3. [Data Structure](#data-structure)
4. [Step-by-Step Analysis](#analysis)
5. [Figure Interpretation (with Examples)](#figures)
6. [Complete Code Examples](#code)
7. [Learning Exercises](#exercises)

## 1. Introduction

### What is This Tutorial?

This tutorial teaches you how to analyze**multi-omics data**from a liver disease study. You'll learn to:
- ✓ Load and explore omics datasets (metabolites, lipids, proteins)
- ✓ Run statistical analyses (linear mixed models)
- ✓ Create 7 types of publication-quality figures
- ✓ Interpret biological results
- ✓ Understand network biology

### The Study

Based on the**MicrobLiver alcohol challenge study**(Israelsen et al., JHEP Reports 2021):

#### Study Parameters:
-**39 participants**: 14 ALD, 15 NAFLD, 10 Healthy controls
-**Alcohol challenge**: Administered via nasogastric tube over 30 minutes
-**3 time points**: 0 (baseline), 60, and 180 minutes after alcohol
-**2 sampling sites**: Hepatic (liver blood) and Systemic (peripheral blood)
-**Multi-omics**: 34 metabolites, 283 lipids, 227 proteins
-**Total observations**: 39 × 3 × 2 = 234 measurements

## 2. Study Design

### Visual Overview

```

39 Participants
├── 14 ALD (Alcohol-related Liver Disease)
├── 15 NAFLD (Metabolic-Associated Fatty Liver Disease)
└── 10 Healthy controls

        ↓ Alcohol Challenge (nasogastric tube, 30 min)

Blood sampling at 3 time points:
├── T0:  0 min (Baseline)
├── T1: 60 min (Early response)
└── T2: 180 min (Late response)

From 2 sites per time point:
├── Hepatic (liver vein) - liver-specific metabolism
└── Systemic (peripheral vein) - whole-body effects

Total: 39 patients × 3 times × 2 sites = 234 observations

```
-**Longitudinal**: Following the same people over time controls for individual differences
-**Dual sites**: Can identify liver-specific vs. systemic changes
-**Multi-omics**: Different molecular layers reveal different aspects of biology

## 3. Data Structure

### File Format

Each dataset is a**data frame**with 234 rows (observations) and the following columns:

| Column | Type | Description | Values |
| --- | --- | --- | --- |
| Specimenno | Character | Unique sample identifier | "ALD_001_hepatic_T0" |
| SampleID | Character | Sample ID | "S1_hepatic_0" |
| PatientID | Character | Patient identifier | "ALD_001", "NAFLD_005", "Healthy_003" |
| pheno | Factor | Disease phenotype | "Healthy", "ALD", "NAFLD" |
| time | Factor | Time point (minutes) | "0", "60", "180" |
| sampleSite | Factor | Sampling location | "hepatic", "systemic" |
| split_col | NA | Separator column | All NA (marks end of metadata) |
| Feature columns | Numeric | Molecular measurements | Scaled concentrations (log10 transformed) |

## 4. Step-by-Step Analysis

### Step 1: Load and Validate Data

```
# Load the data
data <- readRDS("preprocessData/plasma_metabolites_matrix.rds")

# Check structure
str(data)
dim(data)  # Should be 234 x 41

# Validate factor levels
levels(data$pheno)      # "Healthy" "ALD" "NAFLD"
levels(data$time)       # "0" "60" "180"
levels(data$sampleSite) # "hepatic" "systemic"

# View first few rows
head(data[, 1:10])
```

### Step 2: Run Statistical Analysis

```
# Load functions
source("functions.R")

# Set directories
data_dir <- "preprocessData/"
results_dir <- "results/"
models_dir <- "models/"

# Create directories
dir.create(results_dir, showWarnings = FALSE)
dir.create(models_dir, showWarnings = FALSE)

# Run Linear Mixed Models
results <- LMM_time(
  name = "plasma_metabolites",
  ID = "PatientID",
  fixedEffect = "time",
  groupBy = "pheno",
  data_dir = data_dir,
  results_dir = results_dir,
  models_dir = models_dir,
  sampleSite = "hepatic",
  overwrite = TRUE
)

# Check results
head(results)
cat("Significant features (padj<0.05):", sum(results$padj < 0.05, na.rm=TRUE))
```

For each metabolite, fits a model: `metabolite ~ time + (1|PatientID)`
-**Fixed effect (time)**: Tests if metabolite changes at 60 and 180 min vs baseline
-**Random effect (PatientID)**: Accounts for repeated measurements from same patient
-**Grouped by pheno**: Separate analysis for Healthy, ALD, and NAFLD

## 5. Figure Interpretation with Examples

### Figure 1: PCA Plot - Group Separation

![PCA Plot](example_figures/figure1_pca.png)

Figure 1: Principal Component Analysis of Metabolites
-**Each point**= one sample (patient at specific time)
-**Colors**: Green = Healthy, Red = ALD, Orange = NAFLD
-**Shapes**: Circle = 0 min, Square = 60 min, Triangle = 180 min
-**PC1 (x-axis)**: First principal component (explains most variance)
-**PC2 (y-axis)**: Second principal component
- ✓**Group separation**: Do colors cluster separately?
- ✓**Time progression**: Do shapes move in a direction?
- ✓**Outliers**: Any points far from their group?

### Figure 2: Volcano Plot - Statistical Results

![Volcano Plot](example_figures/figure2_volcano.png)

Figure 2: Volcano Plot - ALD vs Healthy at 180 min
-**X-axis**: Effect size (fold change) - how much the metabolite changes
-**Y-axis**: -log10(p-value) - statistical significance (higher = more significant)
-**Red points**: Significantly upregulated (increased in ALD)
-**Blue points**: Significantly downregulated (decreased in ALD)
-**Gray points**: Not significant (p ≥ 0.05)
-**Horizontal line**: Significance threshold (p = 0.05)
-**Vertical line**: No change (fold change = 0)
- ✓**Top-right quadrant**: Significantly increased metabolites
- ✓**Top-left quadrant**: Significantly decreased metabolites
- ✓**Labeled points**: Most significant features
- ✓**Number of significant points**: How many metabolites change?

### Figure 3: Heatmap - Pattern Clustering

![Heatmap](example_figures/figure3_heatmap.png)

Figure 3: Heatmap of Top 20 Significant Metabolites
-**Rows**: Individual metabolites (top 20 most significant)
-**Columns**: Conditions (phenotype × time combinations)
-**Colors**: Red = high concentration, Blue = low concentration, White = medium
-**Dendrograms**: Tree diagrams showing clustering
- ✓**Row clusters**: Metabolites with similar patterns (may be in same pathway)
- ✓**Column clusters**: Conditions with similar metabolic profiles
- ✓**Color blocks**: Coordinated changes across multiple metabolites

### Figure 4: Time Course - Temporal Dynamics

![Time Course](example_figures/figure4_timecourse.png)

Figure 4: Time Course of Top Variable Metabolites
-**X-axis**: Time in minutes (0, 60, 180)
-**Y-axis**: Metabolite concentration (scaled)
-**Lines**: Mean trajectory for each phenotype
-**Shaded areas**: Standard error (confidence interval)
-**Colors**: Green = Healthy, Red = ALD, Orange = NAFLD
- ✓**Rising lines**: Metabolite increases after alcohol
- ✓**Declining lines**: Metabolite decreases after alcohol
- ✓**Diverging lines**: Different responses between phenotypes
- ✓**Parallel lines**: Similar temporal dynamics across groups

### Figure 5: Boxplot - Distribution Analysis

![Boxplot](example_figures/figure5_boxplot.png)

Figure 5: Boxplots of Significant Metabolites
-**Boxes**: Show quartiles (25th, 50th/median, 75th percentile)
-**Whiskers**: Extend to minimum and maximum values
-**Individual dots**: Each data point
-**X-axis**: Time points (0, 60, 180 minutes)
-**Colors**: Different phenotypes
- ✓**Box position**: Median level of metabolite
- ✓**Box height**: Variability within group
- ✓**Overlapping boxes**: Similar distributions
- ✓**Separated boxes**: Significant differences

### Figure 6: Network Graph - Molecular Relationships

![Network by Phenotype](example_figures/figure6_network_by_phenotype.png)

Figure 6: Metabolite Correlation Networks by Phenotype at Baseline
-**Nodes (circles)**: Individual metabolites
-**Node size**: Number of connections (degree) - larger = more connected
-**Node color**: Phenotype (Green=Healthy, Red=ALD, Orange=NAFLD)
-**Edges (lines)**: Significant correlations between metabolites
-**Edge color**: Red = positive correlation, Blue = negative correlation
-**Edge thickness**: Correlation strength (thicker = stronger)
- ✓**Network density**: Many edges = coordinated metabolism
- ✓**Hub nodes**: Large nodes = key regulatory molecules
- ✓**Modules**: Clusters of connected nodes = functional groups
- ✓**Phenotype differences**: Compare network structures
-**Healthy**: 7 nodes, 4 edges - Simple, coordinated network
-**ALD**: 2 nodes, 1 edge -**MOST DISRUPTED**(fewest connections!)
-**NAFLD**: 13 nodes, 9 edges - Most complex network**Interpretation:**ALD shows the most severe metabolic disruption with loss of molecular coordination. NAFLD maintains complexity but with altered connections. This suggests different disease mechanisms: ALD causes metabolic breakdown, while NAFLD causes metabolic rewiring.

### Figure 7: Bar Chart - Summary Statistics

![Bar Chart](example_figures/figure7_barchart.png)

Figure 7: Number of Significant Metabolites at 180 min
-**X-axis**: Comparison groups (ALD vs Healthy, NAFLD vs Healthy)
-**Y-axis**: Number of significant metabolites
-**Red bars**: Upregulated metabolites (increased)
-**Blue bars**: Downregulated metabolites (decreased)
-**Numbers on bars**: Exact counts
- ✓**Bar heights**: Magnitude of metabolic changes
- ✓**Red vs blue balance**: Direction of changes
- ✓**Between-group comparison**: Which disease has more changes?

## 6. Complete Code Example

### Full Analysis Script

```
# ============================================================
# COMPLETE MULTI-OMICS ANALYSIS SCRIPT
# ============================================================

# 1. SETUP
library(lmerTest)
library(caret)
library(ggplot2)
library(ggrepel)

source("functions.R")
source("network_functions.R")

# Directories
data_dir <- "preprocessData/"
results_dir <- "results/"
models_dir <- "models/"
figures_dir <- "figures/"

# Create directories
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(models_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# 2. LOAD DATA
cat("\n=== LOADING DATA ===\n")
metabolites <- readRDS(paste0(data_dir, "plasma_metabolites_matrix.rds"))
cat("Loaded:", nrow(metabolites), "samples\n")

# 3. RUN ANALYSIS
cat("\n=== RUNNING LINEAR MIXED MODELS ===\n")
results <- LMM_time(
  name = "plasma_metabolites",
  ID = "PatientID",
  fixedEffect = "time",
  groupBy = "pheno",
  data_dir = data_dir,
  results_dir = results_dir,
  models_dir = models_dir,
  sampleSite = "hepatic",
  overwrite = TRUE
)

n_sig <- sum(results$padj < 0.05, na.rm=TRUE)
cat("Significant features:", n_sig, "\n")

# 4. CREATE VOLCANO PLOTS
cat("\n=== CREATING VOLCANO PLOTS ===\n")
plots <- volcanoPlot_time(results, pvalue = "padj")

for(name in names(plots)) {
  ggsave(paste0(figures_dir, "volcano_", name, ".png"),
         plots[[name]], width = 8, height = 6, dpi = 300)
}

# 5. CREATE NETWORKS
cat("\n=== CREATING NETWORKS ===\n")
matrix_list <- list(
  plasma_metabolites = readRDS(paste0(data_dir, "plasma_metabolites_matrix.rds")),
  plasma_lipidomics = readRDS(paste0(data_dir, "plasma_lipidomics_matrix.rds")),
  plasma_proteomics = readRDS(paste0(data_dir, "plasma_proteomics_matrix.rds"))
)

for(pheno in c("Healthy", "ALD", "NAFLD")) {
  network <- getMultiOmicsNetwork(
    matrix_list, pheno = pheno, time = "0",
    cor_cutoff = 0.7, figures_path = figures_dir, save = TRUE
  )
}

cat("\n=== COMPLETE! ===\n")
```

## 7. Learning Exercises

### Exercise 1: Compare Sampling Sites

```
# Run analysis on systemic samples
results_systemic <- LMM_time(
  name = "plasma_metabolites",
  ID = "PatientID",
  fixedEffect = "time",
  groupBy = "pheno",
  data_dir = data_dir,
  results_dir = results_dir,
  models_dir = models_dir,
  sampleSite = "systemic",
  overwrite = TRUE
)

# Compare
n_sig_hepatic <- sum(results_hepatic$padj < 0.05, na.rm=TRUE)
n_sig_systemic <- sum(results_systemic$padj < 0.05, na.rm=TRUE)

cat("Hepatic:", n_sig_hepatic, "\n")
cat("Systemic:", n_sig_systemic, "\n")
```**Question:**Do you find more changes in hepatic or systemic samples? What does this tell you about liver-specific metabolism?

### Exercise 2: Analyze Lipids

```
# Run same analysis on lipidomics data
results_lipids <- LMM_time(
  name = "plasma_lipidomics",
  ID = "PatientID",
  fixedEffect = "time",
  groupBy = "pheno",
  data_dir = data_dir,
  results_dir = results_dir,
  models_dir = models_dir,
  sampleSite = "hepatic",
  overwrite = TRUE
)

# Compare with metabolites
cat("Significant metabolites:", sum(results_metabolites$padj < 0.05, na.rm=TRUE), "\n")
cat("Significant lipids:", sum(results_lipids$padj < 0.05, na.rm=TRUE), "\n")
```**Question:**Do lipids or metabolites show more changes? Why might this be?

## 8. Student Checklist
- I can load the synthetic datasets in R
- I understand the study design (39 patients, 3 times, 2 sites)
- I can run the validation script successfully
- I can execute the LMM_time analysis
- I can create volcano plots
- I can create network graphs
- I can interpret a volcano plot
- I can interpret a network graph
- I understand what a linear mixed model is
- I understand why multiple testing correction is needed
- I can identify significant results (padj < 0.05)
- I can compare hepatic vs systemic results
- I can compare different phenotypes
- I can explain the biological meaning of results

## 9. Key Takeaways

### 🎯 Main Findings from This Analysis:
1.**Metabolic Changes**: Both ALD and NAFLD show significant metabolic alterations after alcohol challenge
2.**Network Disruption**: ALD has the most disrupted metabolic network (only 2 nodes, 1 edge)
3.**NAFLD Complexity**: NAFLD maintains network complexity but with altered connections
4.**Temporal Dynamics**: Changes progress over 180 minutes, with different trajectories per phenotype
5.**Disease Mechanisms**: Different network patterns suggest distinct pathological mechanisms

## 10. References**Original Study:**Israelsen M, Kim M, Suvitaival T, et al. Comprehensive lipidomics reveals phenotypic differences in hepatic lipid turnover in ALD and NAFLD during alcohol intoxication. *JHEP Reports*. 2021;3(5):100325.

## 11. Next Steps
- Run the analysis on all omics layers (metabolites, lipids, proteins)
- Create your own custom plots
- Try different correlation thresholds in networks
- Compare networks at different time points
- Explore pathway enrichment analysis

## 🎉 Congratulations!

You now have the skills to analyze multi-omics longitudinal data!

*Keep practicing with different datasets and parameters*