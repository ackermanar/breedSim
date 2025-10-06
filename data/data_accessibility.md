# Data Directory

This directory contains the data structure for the Experimental Design in Plant Breeding project, but **NO DATA FILES ARE STORED IN THIS REPOSITORY - TO BE MADE PUBLIC AFTER MANUSCRIPT RELEASE**. 

## Download Required Data

All data files must be downloaded from cloud storage using the following commands:

### Option 1: Download Individual Files

For specific files, use wget or curl:

**Using wget:**
```bash
wget https://[your-bucket-url]/experimental_design/filename
```

**Using curl:**
```bash
curl -O https://[your-bucket-url]/experimental_design/filename
```

Replace `filename` with the specific file you need from the directory structure below.

### Option 2: Download Entire Dataset via AWS CLI (Recommended)

For the complete dataset, use AWS CLI:

```bash
aws s3 sync s3://[your-bucket]/experimental_design/ ./experimental_design
```

**AWS CLI Setup Required:** Install and configure AWS CLI first. See instructions at:
https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-welcome.html

*Note: No AWS credentials needed for this public dataset - use `aws configure` with dummy values if prompted.*

### Option 3: Download Entire Dataset via Google Cloud CLI

Alternatively, use Google Cloud Storage with gsutil:

```bash
gsutil -m cp -r \
  "gs://[your-bucket]/experimental_design" \
  .
```

**Google Cloud CLI Setup Required:** Install and configure gcloud CLI first. See instructions at:
https://cloud.google.com/storage/docs/discover-object-storage-gcloud

*Note: No authentication required for this public dataset.*

## Directory Structure

After downloading, your data directory should look like this:

```
data/
├── raw/
│   ├── genotypic/
│   │   ├── geno.RData                    # Genotype matrix (4,102 lines × 9,262 SNPs)
│   │   └── K2.RData                      # Genomic relationship matrix
│   └── simulation_parameters/
│       ├── replication_levels.csv        # Replication schemes for trial designs
│       └── cohort_structure.csv          # Breeding cohort definitions (S1-S4)
│
├── simulation_outputs/
│   ├── BLUP/
│   │   ├── results/                      # Individual iteration results (450 runs)
│   │   ├── errors/                       # Error logs from failed simulations
│   │   └── final/                        # Consolidated BLUP results
│   │       └── breedSimResult_BLUP_1-16.csv
│   │
│   ├── gBLUP/
│   │   ├── results/                      # Individual iteration results (50 runs)
│   │   ├── errors/                       # Error logs from failed simulations
│   │   └── final/                        # Consolidated gBLUP results
│   │       └── breedSimResult_gBLUP_allIterationsV2_12-21.csv
│   │
│   └── GEBV/
│       ├── results/                      # Individual iteration results (100 runs)
│       ├── errors/                       # Error logs from failed simulations
│       └── final/                        # Consolidated GEBV sparse testing results
│           └── breedSimResult_sparse_1-16.csv
│
├── processed/
│   ├── Tables/
│   │   ├── Table1_cohort_descriptions.csv
│   │   ├── Table2_trial_designs.csv
│   │   ├── Table3_ANOVA_conventional.csv
│   │   ├── Table4_ANOVA_genomic.csv
│   │   └── Table5_ANOVA_sparse_testing.csv
│   │
│   ├── Figures/
│   │   ├── Fig1_simulation_framework/     # Conceptual diagram data
│   │   ├── Fig2_PCA_population/          # Population structure data
│   │   ├── Fig3_conventional_accuracy/    # BLUP accuracy across parameters
│   │   ├── Fig4_DiD_conventional/        # Difference-in-differences analysis
│   │   ├── Fig5_genomic_accuracy/        # GBLUP accuracy data
│   │   ├── Fig6_DiD_genomic/             # DiD for genomic selection
│   │   ├── Fig7_sparse_testing/          # GEBV accuracy in sparse testing
│   │   └── Fig8_DiD_sparse/              # DiD for sparse testing scenarios
│   │
│   └── Supplemental/
│       ├── SupplTable1_marker_summary.csv
│       ├── SupplTable2_genetic_correlations.csv
│       ├── SupplFig1_convergence_diagnostics.pdf
│       └── SupplFig2_residual_diagnostics.pdf
│
└── scripts/
    ├── breedSimV7_BLUP_1Phase_1-16.R     # BLUP simulation pipeline
    ├── breedSimV9_gBLUP_1-26.R           # gBLUP simulation pipeline
    └── breedSimV9_GEBV_1-26.R            # GEBV sparse testing pipeline
```

## Dataset Details

### Raw Genotypic Data
- **geno.RData**: Genotype-by-sequencing (GBS) marker data for 4,102 wheat breeding lines from University of Illinois breeding program
  - 9,262 SNP markers (MAF > 0.05, call rate > 80%)
  - Lines from four consecutive breeding cohorts (2019-2022)
  - Stages: S1 (preliminary, year 1), S2 (preliminary, year 2), S3 (advanced, year 3), S4 (advanced, year 4)

- **K2.RData**: Genomic relationship matrix calculated using method of Endelman (2011)

### Simulation Parameters
The study simulated multi-environment trials under varying conditions:
- **Heritability (h²)**: 0.2, 0.4, 0.6, 0.8
- **Inter-environment genetic correlation (rGE)**: Fixed at 0.5
- **Intra-environment genetic correlation (localized rGE)**: 0.2, 0.4, 0.6, 0.8, 1.0
- **Replication levels**: Low (1-1-2-2), Intermediate (2-2-2-2), High (2-2-3-3 for S1-S2-S3-S4)
- **Number of environments**: 5
- **Design types**: 
  - Restricted Randomization (RR): Separate preliminary and advanced trials
  - Complete Randomization (CR): All cohorts in single trial using p-rep designs

### Simulation Outputs

#### BLUP Analysis (450 total iterations)
- 500 simulation runs per scenario for conventional selection
- 120 unique parameter combinations (h², rGE, replication level)
- Results include selection accuracy (r_bv) by cohort, trial, and overall

#### gBLUP Analysis (50 total iterations)
- 50 simulation runs per scenario for genomic selection
- Low replication scenarios only
- 40 unique parameter combinations (h², rGE fixed at low replication)
- Includes genomic relationship matrix in mixed models

#### GEBV Sparse Testing (100 total iterations)  
- 50 simulation runs per scenario
- Preliminary cohorts (S1, S2) tested in only 1 of 5 environments
- Advanced cohorts (S3, S4) tested in all 5 environments
- 20 unique parameter combinations (h², rGE)

## Essential Files for Basic Analysis

Download these key files first to replicate the main analyses:

1. `raw/genotypic/geno.RData` - Genotype matrix for simulations
2. `raw/genotypic/K2.RData` - Genomic relationship matrix for gBLUP/GEBV
3. `simulation_outputs/BLUP/final/breedSimResult_BLUP_1-16.csv` - Conventional selection results
4. `simulation_outputs/gBLUP/final/breedSimResult_gBLUP_allIterationsV2_12-21.csv` - Genomic selection results
5. `simulation_outputs/GEBV/final/breedSimResult_GEBV_1-16.csv` - Sparse testing results

## Simulation Framework

The simulations were conducted using ASReml-R 4.2 with the following workflow:

1. **Generate correlated marker effects** across environments and areas (blocks/trials) using Kronecker product of genetic covariance matrices
2. **Calculate breeding values** for each line in each environment-area combination
3. **Simulate phenotypes** by adding environment effects, block effects, and residual error scaled to target heritability
4. **Fit mixed models** with environment-specific variances for block and genetic effects
5. **Estimate breeding values** using BLUP (pedigree only), GBLUP (with genomic relationships), or GEBV (sparse testing)
6. **Calculate selection accuracy** as Pearson correlation between estimated and true breeding values

## Statistical Models

### BLUP (Conventional Selection)
```
y_ijk = μ + e_i + b_k(i) + g_j(i) + ε_ijk
```
Where genetic effects assumed independent (identity matrix)

### GBLUP (Genomic Selection)
```
y_ijk = μ + e_i + b_k(i) + g_j(i) + ε_ijk
```
Where genetic effects follow G_ge = G_e ⊗ G_m with G_m as genomic relationship matrix

### GEBV (Sparse Testing)
Same as GBLUP but with preliminary cohorts phenotyped in only 1 of 5 environments

## Difference-in-Differences (DiD) Analysis

DiD estimation was used to quantify differential responses between RR and CR as conditions deteriorated from optimal:

```
E(Y | Z, T) = β₀ + β₁I(T=1) + β₂Z + β₃I(T=1)Z
```

Where:
- Y = selection accuracy (r_bv)
- Z = design indicator (0=RR control, 1=CR intervention)
- T = condition indicator (0=optimal baseline, 1=suboptimal)
- β₃ = DiD estimator (differential response)

Optimal baseline: h²=0.8, high replication, localized rGE=1.0

## File Naming Conventions

- Iteration results: `results_[iteration]_h[heritability]_[date].csv`
- Error logs: `errors_[iteration]_h[heritability]_[date].Rdata`
- Final consolidated: `breedSimResult_[model]_[date_range].csv`

## Column Descriptions

### Result Files
- **group**: Level of aggregation (overall, cohort, trial, testSet)
- **cor**: Pearson correlation between estimated and true breeding values
- **model**: Analysis method (BLUP, gBLUP, GEBV, GEBVNoPrelim)
- **design**: Randomization scheme (RCBD, PREP)
- **heritability**: Simulated narrow-sense heritability (0.2-0.8)
- **nLoc**: Number of environments (always 5)
- **macroGxE**: Inter-environment genetic correlation (always 0.5)
- **microGxE**: Intra-environment genetic correlation (0.2-1.0)
- **repCat**: Replication category (1-1-2-2, 2-2-2-2, 2-2-3-3)
- **iteration**: Simulation run number

## Computational Requirements

### Software Dependencies
- R ≥ 4.0.0
- ASReml-R 4.2 (commercial license required)
- Required R packages: data.table, tidyverse, purrr, furrr, future, caret, MASS, matrixcalc, MBESS

### Hardware Specifications
- BLUP simulations: ~6 cores, 4GB RAM per worker
- gBLUP/GEBV simulations: ~4-6 cores, 8GB RAM per worker
- Total computational time: ~500 hours across all simulations

## Questions?

For data access issues or questions about the simulation framework:
- **Contact**: Arlyn Ackerman (aja294@cornell.edu) or Jessica Rutkoski (jrutkoski@illinois.edu)
- **Manuscript**: "Randomization across breeding cohorts improves accuracy of conventional and genomic selection"
- **Institution**: University of Illinois at Urbana-Champaign, Department of Crop Sciences

## Citation

If using this dataset, please cite:

Ackerman, A.J. & Rutkoski, J. (2025). Randomization across breeding cohorts improves accuracy of conventional and genomic selection. *[Journal Name]*, *[Volume]*(*Issue*), [pages]. https://doi.org/[DOI]

⚠️ **Important**: The R scripts expect data files to be in these exact directory locations. Ensure files are downloaded to the correct subdirectories before running analyses.
```

---

**Key Adaptations Made:**

1. **Project-Specific Structure**: Organized around BLUP/gBLUP/GEBV simulation outputs rather than metabolomics data
2. **Breeding-Specific Details**: Included cohort structures, replication schemes, and experimental design parameters
3. **Simulation Framework**: Added comprehensive description of the simulation methodology and statistical models
4. **File Organization**: Structured around 450 BLUP, 50 gBLUP, and 100 GEBV simulation iterations
5. **Technical Requirements**: Included computational specifications and software dependencies
6. **Statistical Models**: Provided equations for BLUP, GBLUP, GEBV, and DiD analyses
7. **Column Descriptions**: Defined key variables in the output datasets

This README provides researchers with everything needed to download, understand, and replicate your experimental design simulation study. The structure mirrors the SAP Metabolomics example while being fully adapted to your plant breeding experimental design research.
