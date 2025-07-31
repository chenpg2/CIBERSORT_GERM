# Enhanced Cell Type Deconvolution Analysis Script

## Overview

`deconvolution_germcell.r` is an enhanced germ cell type deconvolution analysis script based on the CIBERSORT algorithm, specifically designed for cell type proportion estimation in testicular tissues.

## Key Features

### 1. Algorithm Improvements
- **CIBERSORT-like Algorithm**: Uses Support Vector Regression (SVR) for cell proportion estimation
- **Feature Gene Selection**: Automatically selects the most discriminative genes for deconvolution
- **TPM Normalization**: Performs TPM normalization and log2 transformation on raw count data
- **Statistical Significance Testing**: Evaluates statistical significance of results through permutation testing

### 2. Cell Type Definitions
The script includes detailed marker genes for 11 cell types:

#### Somatic Cell Types
- **Leydig Cells**: Interstitial cells responsible for testosterone secretion
- **Sertoli Cells**: Supporting cells that maintain the spermatogenesis microenvironment
- **Peritubular Cells**: Peritubular myoid cells

#### Germ Cell Subtypes (Focus)
- **SSC**: Spermatogonial stem cells
- **Diff_Spg**: Differentiating spermatogonia
- **Spermatocytes**: Spermatocytes
- **Round_Spermatids**: Round spermatids
- **Elongating_Spermatids**: Elongating spermatids

#### Immune Cell Types
- **Macrophages**: Macrophages
- **M1_Macrophages**: M1-type macrophages
- **M2_Macrophages**: M2-type macrophages

## Usage

### Basic Usage
```bash
Rscript deconvolution_germcell.r -i input_counts.csv -o output_proportions.csv
```

### Full Parameters
```bash
Rscript deconvolution_germcell.r \
  -i input_counts.csv \
  -o output_proportions.csv \
  -n 0.5 \
  -p 100
```

### Test Mode
```bash
Rscript deconvolution_germcell.r --test
```

## Parameter Description

| Parameter | Long Parameter | Type | Default | Description |
|-----------|----------------|------|---------|-------------|
| -i | --input | String | Required | Input count matrix file (CSV format) |
| -o | --output | String | celltype_proportions.csv | Output cell type proportions file |
| -n | --nu | Numeric | 0.5 | Nu parameter for SVR (0-1) |
| -p | --perm | Integer | 100 | Number of permutations for significance testing |
| -t | --test | Boolean | FALSE | Run with test data |

## Input File Format

The input file should be in CSV format, containing:
- First column: Gene names (as row names)
- Remaining columns: Sample expression data
- Supports raw count data (will be automatically TPM normalized)

Example format:
```
,Sample1,Sample2,Sample3
Gene1,100,150,120
Gene2,50,75,60
...
```

## Output Files

### 1. Main Results File (`*_proportions.csv`)
Contains estimated proportions of each cell type in each sample:
```
Sample,Leydig,Sertoli,SSC,Diff_Spg,...
Sample1,0.15,0.20,0.10,0.08,...
Sample2,0.12,0.18,0.12,0.09,...
```

### 2. Detailed Results File (`*_detailed.csv`)
Contains proportions, p-values, and quality metrics:
- Cell type proportions
- Statistical significance p-values (`*_pval` columns)
- RMSE quality metrics

## Algorithm Principles

### 1. Data Preprocessing
- TPM normalization: `TPM = (counts / gene_length) / sum(RPK) * 1e6`
- Log2 transformation: `log2(TPM + 1)`
- Gene filtering: Retain genes in the top 90% of expression

### 2. Feature Selection
- Calculate coefficient of variation (CV) for genes
- Calculate F-statistic for between-group differences
- Select top 500 feature genes with highest composite scores

### 3. SVR Deconvolution
- Construct signature matrix (based on marker genes)
- Use linear kernel SVR model to estimate cell proportions
- Apply non-negative constraints and proportion normalization

### 4. Statistical Testing
- Permutation testing to assess significance
- Calculate RMSE as fitting quality metric

## Quality Control

### 1. Output Metrics
- **RMSE**: Root Mean Square Error, smaller values indicate better fitting
- **P-value**: Statistical significance, <0.05 indicates significance
- **Proportion Sum**: Sum of all cell type proportions should equal 1

### 2. Recommended Thresholds
- RMSE < 15: Good fitting
- P-value < 0.05: Statistically significant
- Major cell type proportion > 0.01: Reliable results

## Important Notes

1. **Input Data Quality**: Ensure gene names are standardized, recommend using official gene symbols
2. **Sample Size**: Recommend at least 3 samples for reliable statistical results
3. **Computational Resources**: Permutation testing is time-consuming, can appropriately reduce permutation count
4. **Result Interpretation**: Proportion results are relative quantification, need to be interpreted with biological context

## Dependencies

- `tidyverse`: Data processing
- `optparse`: Command line argument parsing
- `e1071`: Support vector machines
- `parallel`: Parallel computing

## Example Analysis

Run test mode to see example:
```bash
Rscript deconvolution_germcell.r --test
```

This will generate:
- `test_expression_data.csv`: Simulated expression data
- `test_celltype_proportions.csv`: Cell proportion results
- `test_celltype_proportions_detailed.csv`: Detailed results

## Update Log

- **v2.0**: Integrated CIBERSORT algorithm, added statistical testing
- **v1.0**: Basic version, simple averaging method

## Contact Information

Author: chenpeigen  
Date: 2025-07-31
