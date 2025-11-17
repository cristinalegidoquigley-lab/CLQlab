# CLQlab
tutorials from our lab
[INDEX.md](https://github.com/user-attachments/files/23585254/INDEX.md)
# Synthetic MicrobLiver Datasets - Complete Package

## 📦 What You're Getting

I've created synthetic datasets that exactly match your MicrobLiver study structure, allowing you to test and develop your analysis code.

## 📁 Files Included

### Main Datasets (in preprocessData/ folder)
- `plasma_metabolites_matrix.rds` (63 KB) - 34 metabolites
- `plasma_lipidomics_matrix.rds` (501 KB) - 283 lipids  
- `plasma_proteomics_matrix.rds` (139 KB) - 77 proteins
- `plasma_MS_proteomics_matrix.rds` (267 KB) - 150 MS proteins

### Documentation
- `preprocessData/README.md` - Detailed documentation of data structure
- `SUMMARY.md` - Quick overview and usage guide
- `show_data_structure.R` - R script to validate and explore the datasets

## 🔢 Dataset Specifications

**Study Design:**
- 39 participants (14 ALD, 15 NAFLD, 10 Healthy)
- 3 time points (0, 60, 180 minutes)
- 2 sampling sites (hepatic, systemic)
- 234 total observations (39 × 3 × 2)

**Data Properties:**
- All data is preprocessed (log10 transformed and scaled)
- Includes phenotype-specific effects
- Includes time-dependent changes
- Includes sampling site differences
- Random feature names to protect original data

## ✅ Compatibility

These datasets work with ALL functions in your `functions.R` file:
- ✓ `readFile()` - Load and filter data
- ✓ `LMM_time()` - Linear mixed models analysis
- ✓ `volcanoPlot_time()` - Volcano plot visualization
- ✓ `getBoxplot()` - Boxplot generation
- ✓ `getChordDiag()` - Chord diagram for correlations
- ✓ `getCircularHeatmap()` - Circular heatmap visualization

## 🚀 Quick Start

```r
# Load a dataset
data <- readRDS("preprocessData/plasma_metabolites_matrix.rds")

# Check the structure
str(data)
head(data)
table(data$pheno)

# Or run the validation script
source("show_data_structure.R")
```

## 📊 Data Structure

Each dataset has this structure:
```
Columns 1-7: Metadata
  - Specimenno: Sample identifier
  - SampleID: Sample ID
  - PatientID: Patient identifier  
  - pheno: Disease group (Healthy, ALD, NAFLD)
  - time: Time point (0, 60, 180)
  - sampleSite: Location (hepatic, systemic)
  - split_col: Separator (all NA)

Columns 8+: Feature measurements
  - Metabolite_001, Metabolite_002, etc.
  - Lipid_001, Lipid_002, etc.
  - Protein_001, Protein_002, etc.
```

## 💡 Use Cases

Perfect for:
- Testing your analysis pipeline
- Developing new visualization functions
- Training new team members
- Debugging code
- Creating analysis templates
- Preparing for grant applications or presentations

## ⚠️ Important Notes

1. **Not real data**: All values are synthetically generated
2. **Random names**: Feature names are generic placeholders
3. **For testing only**: Not suitable for biological interpretation
4. **Reproducible**: Generated with set.seed(123)

## 📚 Reference

Based on the structure from:
*Israelsen M, et al. Comprehensive lipidomics reveals phenotypic differences in hepatic lipid turnover in ALD and NAFLD during alcohol intoxication. JHEP Reports. 2021;3(5):100325.*

## 🎯 Everything is Ready!

All datasets are validated and ready to use with your existing analysis code. Just load the RDS files and run your analyses!
