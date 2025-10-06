# Randomization Across Breeding Cohorts Improves Accuracy of Conventional and Genomic Selection

[![DOI](https://img.shields.io/badge/DOI-pending-blue)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains simulation code and analysis scripts for evaluating how experimental design impacts prediction accuracy in plant breeding programs. The study demonstrates that randomizing breeding materials across cohorts improves selection accuracy compared to traditional cohort-separated trials, particularly when data is incomplete or sparse. An extensive user guide for user implementation can be found in docs/simulations_parameters_and_specifications_guide.md and .pdf.

**Key Finding**: Complete randomization (CR) across breeding cohorts outperforms restricted randomization (RR) by up to 15.7% when phenotypic or genomic data is limited, though both designs perform equivalently with complete genomic and phenotypic datasets.

## Citation

Ackerman, A.J. & Rutkoski, J. (2025). Randomization across breeding cohorts improves accuracy of conventional and genomic selection. *Manuscript in preparation*.

## Background

### The Problem

Breeding programs typically evaluate materials in separate yield trials based on their advancement stage (cohorts). This spatial separation can confound genetic effects with non-genetic trial effects, potentially reducing selection accuracy—especially when:
- Genomic relationship data is unavailable (conventional BLUP)
- Testing designs are sparse (not all lines in all environments)
- Genotype-by-environment (G×E) interactions are strong

### The Solution

We compared two randomization strategies:
- **Restricted Randomization (RR)**: Traditional approach where cohorts occupy separate trials within environments
- **Complete Randomization (CR)**: Alternative approach where all cohorts are randomized together (e.g., p-rep designs)

## Repository Structure

```
.
├── R/
│   ├── breedSimV7_BLUP_1Phase_1-16.R      # Conventional BLUP simulation
│   ├── breedSimV9_gBLUP_1-26.R            # Genomic BLUP (balanced MET)
│   └── breedSimV9_GEBV_1-26.R             # Genomic-enabled sparse testing
├── data/
│   ├── geno.RData                          # Genotypic data matrix (not included - see Data Availability)
│   └── K2.RData                            # Genomic relationship matrix (generated from geno.RData)
├── docs/
│   ├── simulations_parameters_and_specifications_guide.md            # Detailed parameter specifications
│   └── results_summary.md                                            # Key findings summary
├── output/
│   ├── BLUPoutput/                        # Conventional selection results
│   ├── gBLUPoutput/                       # Genomic selection results (balanced)
│   └── Sparseoutput/                      # Sparse testing results
├── LICENSE
└── README.md
```

## Simulation Framework

### Germplasm

Real wheat breeding lines from the University of Illinois wheat breeding program:
- **4 breeding cohorts** (S1-S4) representing different advancement stages
- **3,102 experimental lines** + check varieties
- **9,262 SNP markers** (GBS, MAF > 0.05, <10% missing data)

### Experimental Design Parameters

**Randomization schemes**:
- Restricted Randomization (RR): RCBD or IBD designs with cohorts in separate trials
- Complete Randomization (CR): p-rep designs with all cohorts randomized together

**Replication levels**:
1. **High**: S4 replicated 3×, S1-S3 replicated 2×
2. **Intermediate**: All cohorts replicated 2×
3. **Low**: S1-S2 unreplicated, S3-S4 replicated 2×

**Environmental parameters**:
- **5 environments** (location × year combinations)
- **Inter-environment genetic correlation**: r = 0.5 (fixed)
- **Intra-environment genetic correlation** (r<sub>GE</sub>): 0.2, 0.4, 0.6, 0.8, 1.0
- **Heritability** (h²): 0.2, 0.4, 0.6, 0.8

### Statistical Models

**Conventional Selection (BLUP)**:
```
y_ijk = μ + e_i + b_k(i) + g_j(i) + ε_ijk
```

**Genomic Selection (GBLUP)**:
```
y_ijk = μ + e_i + b_k(i) + g_j(i) + ε_ijk
where g ~ N(0, G_m ⊗ G_e)
```

**Sparse Testing (GEBV)**:
- S1 and S2 cohorts evaluated in only 1 of 5 environments
- S3 and S4 cohorts evaluated in all 5 environments
- Prediction accuracy assessed for untested cohorts

### Analysis Approach

**Difference-in-Differences (DiD) Framework**:
```
E(Y | Z, T) = β₀ + β₁I(T=1) + β₂Z + β₃I(T=1)Z
```

Where:
- Y = prediction accuracy (r<sub>bv</sub>)
- Z = randomization scheme (0=RR, 1=CR)
- T = parameter conditions (0=optimal, 1=suboptimal)
- β₃ = DiD estimator (differential response between designs)

## Requirements

### Software
- R (≥ 4.0.0)
- ASReml-R (v4.2) - *Requires license*

### R Packages
```r
# Core packages
library(data.table)
library(tidyverse)
library(magrittr)

# Parallel processing
library(furrr)
library(future)

# Statistical modeling
library(asreml4)      # Commercial license required
library(MASS)
library(matrixcalc)
library(MBESS)

# Design generation
library(caret)
library(purrr)
```

## Installation & Usage

### 1. Clone Repository
```bash
git clone https://github.com/yourusername/breeding-cohort-randomization.git
cd breeding-cohort-randomization
```

### 2. Prepare Data
Due to data sharing restrictions, genotypic data is not included. To replicate:

```r
# Load your own marker data (matrix format: lines × markers)
# Rows = genotype names, Columns = SNP markers coded as -1, 0, 1
geno <- your_marker_matrix

# Calculate genomic relationship matrix
library(rrBLUP)
K2 <- A.mat(geno - 1)  # ASReml requires mean-centered markers

# Save for simulation scripts
saveRDS(geno, "data/geno.RData")
saveRDS(K2, "data/K2.RData")
```

### 3. Run Simulations

**Conventional Selection (BLUP)**:
```bash
Rscript R/breedSimV7_BLUP_1Phase_1-16.R
```
- 450 iterations
- Tests all 120 parameter combinations
- Output: `output/BLUPoutput/BLUPresults_final/`

**Genomic Selection (balanced MET)**:
```bash
Rscript R/breedSimV9_gBLUP_1-26.R
```
- 50 iterations
- Low replication scenarios only
- Output: `output/gBLUPoutput/gBLUPresults_final/`

**Genomic Prediction (sparse testing)**:
```bash
Rscript R/breedSimV9_GEBV_1-26.R
```
- 100 iterations per heritability level
- Sparse testing scenarios
- Output: `output/GEBVoutput/GEBVresults_final/`

### 4. Analyze Results

Each simulation outputs CSV files containing:
- `germplasmName`: Genotype identifier
- `cor`: Prediction accuracy (correlation between true and predicted breeding values)
- `design`: RCBD or PREP
- `heritability`: Simulated h²
- `nLoc`: Number of environments
- `macroGxE`: Inter-environment genetic correlation
- `microGxE`: Intra-environment genetic correlation (r<sub>GE</sub>)
- `repCat`: Replication category
- `group`: Overall, by-cohort, or by-test results
- `model`: BLUP, gBLUP, or GEBV
- `iteration`: Simulation replicate number

## Key Results

### 1. Conventional Selection (BLUP)

| Scenario | CR Advantage | Significance |
|----------|--------------|--------------|
| Overall (incomplete data) | +8.3 pp (11.7%) | *** |
| Low replication | +10.1 pp (15.7%) | *** |
| r<sub>GE</sub> = 0.2 | +18.6 pp | *** |

### 2. Genomic Selection (GBLUP - balanced MET)

With complete phenotypic and genomic data:
- **No significant difference** between designs (r<sub>bv</sub> = 0.888 vs 0.885)
- Genomic relationships eliminate confounding effects

### 3. Sparse Testing (GEBV)

| Condition | CR Advantage | DiD coefficient |
|-----------|--------------|-----------------|
| Overall | +1.5% | δ̂ = 0.018 ** |
| h² = 0.2, r<sub>GE</sub> = 0.2 | +5.5% | Highly significant |

### Factors Influencing Design Performance

**Most Critical**: r<sub>GE</sub> (intra-environment genetic correlation)
- DiD: δ̂ = 0.082, p < 0.001 (BLUP)
- Each 0.2 decrease in r<sub>GE</sub> = 8.2 pp advantage for CR

**Moderate Impact**: Replication level
- DiD: δ̂ = 0.005, p < 0.001
- Largest effect when moving to unreplicated entries

**Minimal Impact**: Heritability
- Both designs respond similarly to decreasing h²
- DiD: non-significant across all models

## Practical Recommendations

### Use Complete Randomization (CR) When:
1. **Limited genomic data** or relying on phenotypic BLUP
2. **Implementing sparse testing** designs
3. **Strong G×E interactions** expected (low r<sub>GE</sub>)
4. **Resource constraints** limit replication

### Restricted Randomization (RR) Acceptable When:
1. **Comprehensive genomic data** available for all lines
2. **Balanced multi-environment testing** implemented
3. **Within-cohort selection** is primary objective
4. **Logistical constraints** favor separate trials

## Computational Requirements

| Simulation Type | Iterations | Cores | RAM | Time |
|----------------|-----------|-------|-----|------|
| BLUP | 450 | 6 | 16 GB | ~48 hrs |
| GBLUP | 50 | 4 | 32 GB | ~24 hrs |
| GEBV | 400 | 4 | 32 GB | ~72 hrs |

*Times are approximate and depend on hardware specifications*

## Troubleshooting

### Common Issues

**ASReml convergence failures**:
```r
# Increase workspace
asreml.options(pworkspace = "8gb")

# Adjust AI settings
asreml.options(ai.sing = TRUE, fail = "soft")
```

**Memory issues with large datasets**:
```r
# Use data.table for efficiency
library(data.table)
setDTthreads(threads = 0)  # Use all available threads
```

**Matrix singularity warnings**:
- Ensure sufficient genetic variation in cohorts
- Check for duplicate genotypes
- Verify genomic relationship matrix is positive definite

## Data Availability

Due to data sharing agreements:
- **Genotypic data**: Available upon reasonable request to the authors
- **Simulated phenotypes**: Generated de novo by scripts in this repository
- **Summary statistics**: Included in `docs/results_summary.md`

## Contributing

We welcome contributions! Please:
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/improvement`)
3. Commit changes (`git commit -am 'Add improvement'`)
4. Push to branch (`git push origin feature/improvement`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

## Contact

**Arlyn Ackerman**  
Breeding Insight, Cornell University  
Email: aja258@cornell.edu

**Jessica Rutkoski**  
Department of Crop Sciences, University of Illinois at Urbana-Champaign  
Email: rutkoski@illinois.edu

## Acknowledgments

- Eastern Regional Small Grains Genotyping Lab for GBS services
- Breeding Insight for computational resources

## References

Key methodological references:

1. **Experimental Design**:
   - Cullis et al. (2006) - p-rep designs
   - Clarke & Stefanova (2011) - Optimal designs for early-generation trials
   - Piepho & Williams (2006) - Restricted vs. complete randomization

2. **Genomic Selection**:
   - Combs & Bernardo (2013) - GS accuracy factors
   - Atanda et al. (2022) - Sparse testing with GS

3. **Statistical Methods**:
   - Smith et al. (2007) - Environment-specific variance models
   - Rothbard et al. (2024) - Difference-in-differences methodology

---

**Last Updated**: January 2025  
**Status**: Manuscript in preparation
