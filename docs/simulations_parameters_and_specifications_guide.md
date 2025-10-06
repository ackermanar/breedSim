# Simulation Parameter Specifications Guide

**Version**: 1.1  
**Last Updated**: January 2025  
**Corresponding Manuscript**: "Randomization across breeding cohorts improves accuracy of conventional and genomic selection"

**Version History**:
- v1.1 (2025-01-XX): Corrected r_{k_{1}k_{2}} terminology; changed GEBV to Sparse Testing
- v1.0 (2025-01-XX): Initial comprehensive guide

---

## Table of Contents

1. [Overview](#overview)
2. [Germplasm Structure](#germplasm-structure)
3. [Experimental Design Parameters](#experimental-design-parameters)
4. [Genetic Architecture](#genetic-architecture)
5. [Environmental Parameters](#environmental-parameters)
6. [Phenotype Simulation](#phenotype-simulation)
7. [Statistical Models](#statistical-models)
8. [Simulation Scenarios](#simulation-scenarios)
9. [Output Specifications](#output-specifications)
10. [Parameter Sensitivity](#parameter-sensitivity)

---

## Overview

This guide provides comprehensive documentation of all simulation parameters used to evaluate restricted randomization (RR) versus complete randomization (CR) in multi-environment breeding trials. The simulation framework generates realistic breeding trial data by combining actual wheat genotypic data with simulated phenotypes that reflect quantitative trait architecture, genotype-by-environment interaction, and experimental error.

### Simulation Philosophy

Our approach prioritizes biological realism over computational convenience:
- **Real marker data** from University of Illinois wheat breeding lines
- **Realistic cohort structure** reflecting actual breeding program advancement
- **Biologically plausible variance partitioning** based on empirical wheat yield trial data
- **Conservative parameter ranges** representing typical breeding program conditions

### Simulation Workflow

```
1. Load genotypic data (M matrix: lines × markers)
2. Define cohort structure and trial designs
3. Generate genetic covariance matrices (macro and micro GxE)
4. Simulate marker effects from multivariate normal distribution
5. Calculate breeding values (BV = M × marker effects)
6. Add environmental and spatial effects
7. Add residual error based on target heritability
8. Fit mixed models (BLUP, GBLUP, or Sparse Testing)
9. Calculate prediction accuracy (r_bv)
10. Repeat across parameter combinations and iterations
```

---

## Germplasm Structure

### Cohort Composition

The simulation uses four breeding cohorts representing consecutive advancement stages in a typical wheat breeding program:

| Cohort | Description | N Lines | Designation | Rep Level (typical) |
|--------|-------------|---------|-------------|---------------------|
| S1 | Stage 1, first-year yield testing | 1,204 | Preliminary | Unreplicated |
| S2 | Stage 2, second-year yield testing | 1,519 | Preliminary | 1-2 replicates |
| S3 | Stage 3, third-year yield testing | 1,151 | Advanced | 2 replicates |
| S4 | Stage 4, fourth-year yield testing | 228 | Advanced | 2-3 replicates |
| **Total** | Experimental lines | **4,102** | — | — |

### Check Varieties

In addition to experimental lines, check varieties are included:
- **Function**: Estimate block effects and error variance in designs with unreplicated entries
- **Number**: Variable depending on design (typically 2 per block)
- **Treatment**: Assigned to separate cohort (`cohort = "check"`)

### Genotypic Data Specifications

**Source**: University of Illinois wheat breeding program (2019-2022 cohorts)

**Marker Platform**: Genotyping-by-sequencing (GBS)
- Performed by Eastern Regional Small Grains Genotyping Lab
- Methods described in Gaire et al. (2021, 2022)

**Quality Control Filtering**:
```r
# Initial markers
Raw SNPs: 12,505

# Filtering criteria
Call rate: ≥ 80%
Minor allele frequency (MAF): ≥ 0.05
Missing data: < 10%

# Final marker set
Filtered SNPs: 9,262
```

**Marker Coding**:
- **Homozygous reference**: -1
- **Heterozygous**: 0  
- **Homozygous alternate**: 1
- **Missing data**: Imputed using mean allele frequency

**Genomic Relationship Matrix (GRM)**:
```r
# Calculated using VanRaden's method (Endelman 2011)
K = (M - 1P') (M - 1P')' / (2 Σp_i(1-p_i))

where:
  M = marker matrix (n × m)
  P = matrix of allele frequencies
  p_i = frequency of reference allele at locus i
```

### Population Structure

**Principal Component Analysis Results**:
- **PC1**: 5.8% of variance
- **PC2**: 5.0% of variance
- **Interpretation**: Minimal population structure, as expected for lines from a single breeding program
- **Cohort separation**: Slight differentiation between preliminary (S1, S2) and advanced (S3, S4) stages

**Genetic Relatedness**:
- Lines within cohorts are more closely related than lines between cohorts
- All lines derived from common parent pool used in breeding program
- Average relationship coefficient within cohorts: ~0.10-0.15
- Average relationship coefficient between cohorts: ~0.05-0.10

---

## Experimental Design Parameters

### Design Type Specifications

#### Restricted Randomization (RR)

Cohorts are spatially separated into different trials within each environment.

**Trial Structure**:

| Trial Type | Cohorts | Design | Blocks | Plot Allocation |
|------------|---------|--------|--------|-----------------|
| Preliminary | S1, S2 | RCBD or CRD | 1-2 | ~70% of environment |
| Advanced | S3, S4 | RCBD or IBD | 2-3 | ~30% of environment |

**Key Characteristics**:
- Randomization occurs **within each trial**
- Different trials may have different designs (e.g., RCBD vs IBD)
- Check varieties replicated within each trial separately
- Block effects nested within trial

**Implementation**:
```r
# Preliminary trial
prelimTrial <- breedSim %>%
  filter(cohort %in% c("S1", "S2")) %>%
  mutate(trial = paste0("Loc", location, "_prelim"))

# Advanced trial
advTrial <- breedSim %>%
  filter(cohort %in% c("S3", "S4")) %>%
  mutate(trial = paste0("Loc", location, "_adv"))
```

#### Complete Randomization (CR)

All cohorts are randomized together in a single trial per environment.

**Trial Structure**:

| Trial Type | Cohorts | Design | Blocks | Plot Allocation |
|------------|---------|--------|--------|-----------------|
| Multi-cohort | S1, S2, S3, S4 | p-rep or IBD | 2-3 | 100% of environment |

**Key Characteristics**:
- Randomization occurs **across all cohorts**
- Single unified design per environment
- Replicated lines (S3, S4) serve as "checks" for unreplicated lines (S1, S2)
- Block effects span all cohorts

**Implementation**:
```r
# Single multi-cohort trial
multiCohortTrial <- breedSim %>%
  mutate(trial = paste0("Loc", location, "_PREP"))
```

### Replication Level Specifications

Three replication levels tested to reflect different resource allocation strategies:

#### 1. High Replication

**Rationale**: Programs with dedicated elite trials or abundant resources

| Design | S1 | S2 | S3 | S4 | Checks/Block |
|--------|----|----|----|----|--------------|
| RR - Preliminary | 2 | 2 | — | — | 2 |
| RR - Advanced | — | — | 2 | 3 | 2 |
| CR - Multi-cohort | 2 | 2 | 2 | 3 | 2 |

**Total plots per environment**: ~18,000-20,000

**Statistical power implications**:
- High precision for all cohorts
- Minimal reliance on spatial or genetic borrowing of information
- Suitable for within-cohort selection decisions

#### 2. Intermediate Replication

**Rationale**: Balanced resource allocation typical of many breeding programs

| Design | S1 | S2 | S3 | S4 | Checks/Block |
|--------|----|----|----|----|--------------|
| RR - Preliminary | 2 | 2 | — | — | 2 |
| RR - Advanced | — | — | 2 | 2 | 2 |
| CR - Multi-cohort | 2 | 2 | 2 | 2 | 2 |

**Total plots per environment**: ~16,000-18,000

**Statistical power implications**:
- Adequate precision for all cohorts
- Moderate reliance on spatial models
- Suitable for across-cohort selection and parent selection

#### 3. Low Replication

**Rationale**: Resource-constrained programs or sparse testing strategies

| Design | S1 | S2 | S3 | S4 | Checks/Block |
|--------|----|----|----|----|--------------|
| RR - Preliminary | 1 | 1 | — | — | 2 |
| RR - Advanced | — | — | 2 | 2 | 2 |
| CR - Multi-cohort | 1 | 1 | 2 | 2 | 2 |

**Total plots per environment**: ~10,000-12,000

**Statistical power implications**:
- Unreplicated early-stage entries require spatial or genetic borrowing
- High dependence on replicated entries for error variance estimation
- Requires ≥20% replicated plots for reliable block effect estimation (Clarke & Stefanova 2011)

### Block Structure

#### Restricted Randomization (RR)

Blocks are **nested within trials**:

```r
# Block assignment for RR
prelimBlocks <- max(breedSim[cohort %in% c("S1", "S2"), repNum])
advBlocks <- max(breedSim[cohort %in% c("S3", "S4"), repNum])

breedSim[test == "prelim", block := repNum]
breedSim[test == "adv", block := repNum + prelimBlocks]
```

**Example Structure (Intermediate Replication)**:
```
Environment 1
├── Preliminary Trial
│   ├── Block 1 (S1, S2, checks - rep 1)
│   └── Block 2 (S1, S2, checks - rep 2)
└── Advanced Trial
    ├── Block 3 (S3, S4, checks - rep 1)
    └── Block 4 (S3, S4, checks - rep 2)
```

#### Complete Randomization (CR)

Blocks are **resolvable and span all cohorts**:

```r
# Block assignment for CR
totalBlocks <- max(breedSim[cohort == "S4", repTotal])

# Fully replicated cohorts assigned directly
breedSim[cohort == "S4", block := repNum]
breedSim[cohort == "S3", block := repNum]

# Partially replicated cohorts distributed across blocks
if (S2reps < totalBlocks) {
  folds <- createFolds(breedSim[cohort == "S2", germplasmName], 
                       k = totalBlocks, list = FALSE)
  breedSim[cohort == "S2", block := folds]
}
```

**Example Structure (Intermediate Replication)**:
```
Environment 1 - Single Multi-cohort Trial
├── Block 1 (S1_subset, S2_subset, S3_rep1, S4_rep1, checks)
├── Block 2 (S1_subset, S2_subset, S3_rep2, S4_rep2, checks)
```

### Check Plot Allocation

**Restricted Randomization**:
- Checks replicated within each trial independently
- Preliminary trial: 2 × preliminary blocks
- Advanced trial: 2 × advanced blocks
- Total checks = 2 × (preliminary blocks + advanced blocks)

**Complete Randomization**:
- Checks distributed across all blocks
- Total checks = 2 × total blocks
- More efficient allocation (fewer total check plots)

**Example Calculation (5 environments, intermediate replication)**:

| Design | Prelim Blocks | Adv Blocks | Total Blocks | Check Plots/Env |
|--------|---------------|------------|--------------|-----------------|
| RR | 2 | 2 | 4 | 8 (2×2 + 2×2) |
| CR | — | — | 2 | 4 (2×2) |

---

## Genetic Architecture

### Marker Effect Simulation

Marker effects are simulated to create realistic quantitative trait architecture with genotype-by-environment interaction at two spatial scales.

#### Genetic Covariance Structure

**Two-level GxE model**:

1. **Macro-GxE**: Genetic correlation between environments (location × year)
2. **Intra-environment GxE**: Genetic correlation between areas within environments (r_{k_{1}k_{2}})

**Covariance Matrix Construction**:

The genetic covariance matrix is constructed using the Kronecker product:

```
Σ_g = J ⊗ K

where:
  J = macro-GxE matrix (environment × environment)
  K = intra-environment GxE matrix (area × area within environment)
  ⊗ = Kronecker product
```

**Matrix J (Macro-GxE)**:

Fixed at inter-environment correlation of **r = 0.5** across all simulations.

For 5 environments:
```
J = [1.0  0.5  0.5  0.5  0.5]
    [0.5  1.0  0.5  0.5  0.5]
    [0.5  0.5  1.0  0.5  0.5]
    [0.5  0.5  0.5  1.0  0.5]
    [0.5  0.5  0.5  0.5  1.0]
```

**Rationale for r = 0.5**:
- Based on empirical wheat yield trial correlations (Lado et al. 2016)
- Represents moderate G×E typical of regional testing
- Allows genetic signal to transfer across environments while maintaining realistic heterogeneity

**Matrix K (Intra-environment GxE)**:

**Varied as treatment parameter**: r_{k_{1}k_{2}} = 0.2, 0.4, 0.6, 0.8, 1.0

For RR with 2 trials per environment (preliminary, advanced):
```
K_RR = [1.0           r_{k_{1}k_{2}}]
       [r_{k_{1}k_{2}}   1.0        ]
```

For CR with 2 blocks per environment:
```
K_CR = [1.0           r_{k_{1}k_{2}}]
       [r_{k_{1}k_{2}}   1.0        ]
```

**Interpretation of r_{k_{1}k_{2}} values**:

| r_{k_{1}k_{2}} | Biological Interpretation | Breeding Implications |
|----------------|---------------------------|----------------------|
| 1.0 | No area-specific GxE | Ranks identical across trials/blocks |
| 0.8 | Weak area-specific GxE | Minor rank changes |
| 0.6 | Moderate area-specific GxE | Noticeable rank changes |
| 0.4 | Strong area-specific GxE | Substantial rank changes |
| 0.2 | Very strong area-specific GxE | Severe rank changes (extreme scenario) |

#### Example: Full Covariance Matrix

For 2 environments, 2 areas per environment, r_{k_{1}k_{2}} = 0.6:

```
Σ_g = J ⊗ K = [1.0  0.5] ⊗ [1.0  0.6]
              [0.5  1.0]   [0.6  1.0]

    = [1.0  0.6  0.5  0.3]  <- Env1-Area1, Env1-Area2, Env2-Area1, Env2-Area2
      [0.6  1.0  0.3  0.5]
      [0.5  0.3  1.0  0.6]
      [0.3  0.5  0.6  1.0]
```

**Element interpretation**:
- Diagonal (1.0): Variance for each area
- Within-environment off-diagonal (0.6): Correlation between areas in same environment
- Between-environment, same area type (0.5): Correlation across environments
- Between-environment, different area type (0.3): Product of macro and intra-environment correlations (0.5 × 0.6)

#### Marker Effect Sampling

**Multivariate normal distribution**:

```r
# Sample marker effects
markerEffects <- mvrnorm(
  n = ncol(geno),              # Number of markers (9,262)
  mu = rep(0, nrow(Σ_g)),      # Mean = 0 for all areas
  Sigma = Σ_g                   # Genetic covariance matrix
)
```

**Output dimensions**: 9,262 markers × (n_environments × n_areas_per_environment)

**For our simulation**:
- 5 environments × 2 areas = 10 area-specific marker effect vectors
- Each marker has 10 different effects representing its performance in each area

#### Positive Definite Matrix Correction

Not all combinations of r_macro and r_{k_{1}k_{2}} produce positive definite matrices. When non-positive definite:

```r
if (!is.positive.definite(Σ_g)) {
  Σ_g <- nearPD(Σ_g, keepDiag = TRUE, corr = TRUE, do2eigen = FALSE)
  Σ_g <- as.matrix(Σ_g[["mat"]])
  warning("Matrix corrected using nearPD")
}
```

**Note**: All parameter combinations used in final simulations produced positive definite matrices without correction.

### Breeding Value Calculation

#### Area-Specific Breeding Values

```r
# Calculate breeding values for each area
BV_area = M × β_area

where:
  M = genotype matrix (n_lines × n_markers)
  β_area = marker effects for specific area (n_markers × 1)
```

**Result**: Each line has different breeding values in different areas, creating realistic G×E interaction.

#### True Breeding Value

The "true" breeding value used for accuracy assessment is the **mean across all areas**:

```r
trueBV_i = mean(BV_i,area1, BV_i,area2, ..., BV_i,area_p)
```

**Rationale**:
- Represents overall genetic merit across target population of environments
- Analogous to breeding for broad adaptation
- Provides single reference value for accuracy calculations

#### Genetic Variance Partitioning

**Within-environment genetic variance**:
```r
σ²_g(env) = var(BV_area)  # Calculated separately for each environment
```

**Between-environment genetic variance**:
```r
σ²_GxE = var(mean(BV_area by environment)) 
```

**Typical variance ratios observed**:
- σ²_g(env) ≈ 1.0 (standardized by sampling process)
- σ²_GxE / σ²_g ≈ 0.3-0.5 (depending on r_{k_{1}k_{2}})

---

## Environmental Parameters

### Number of Environments

**Fixed at n = 5** across all simulations.

**Rationale**:
- Typical of regional yield trial networks
- Sufficient to estimate G×E while remaining computationally feasible
- Represents minimum needed for stable variance component estimation

**Interpretation**: Each environment represents unique location × year combination (e.g., Urbana-2022, DeKalb-2022, etc.)

### Environment Main Effects

Environment main effects represent systematic productivity differences between locations/years.

**Sampling distribution**:
```r
e_j ~ N(0, 2σ²_g(j))

where:
  σ²_g(j) = genetic variance in environment j
```

**Variance ratio**: Environment variance = **2× genetic variance**

**Rationale**:
- Based on empirical wheat yield trials where σ²_env / σ²_g ≈ 2-3
- Ensures environment is major source of phenotypic variation
- Creates realistic partition: genetic < environment effects

**Example values** (h² = 0.6):
```
σ²_g = 1.0
σ²_env = 2.0
σ²_ε = 0.67
σ²_p = 3.67
```

### Area Effects (Block/Trial Effects)

Area effects represent local heterogeneity within environments (soil fertility, moisture gradients, field position).

**Sampling distribution**:
```r
b_k(j) ~ N(0, 0.5σ²_g(j))

where:
  σ²_g(j) = genetic variance in environment j
```

**Variance ratio**: Area variance = **0.5× genetic variance**

**Rationale**:
- Spatial heterogeneity is generally smaller than genetic variance in well-managed trials
- Ratio of 0.5 represents moderate spatial trends
- Comparable to block variance estimates in wheat yield trials

**Design implications**:
- **RR**: Area = trial (preliminary vs advanced)
  - Creates systematic differences between trials beyond just genetic composition
  - Can confound with cohort genetic means
  
- **CR**: Area = resolvable block within unified trial
  - Represents local field heterogeneity
  - No confounding with cohort structure

### Inter-Environment Correlation

**Fixed at r = 0.5** (see Genetic Architecture section)

### Intra-Environment Correlation (r_{k_{1}k_{2}})

**Treatment parameter**: 0.2, 0.4, 0.6, 0.8, 1.0

**Key driver of design performance differences** - see Results section.

---

## Phenotype Simulation

### Phenotype Model

Complete phenotypic model for line *i* in environment *j* and area *k*:

```
y_ijk = g_i(jk) + e_j + b_k(j) + ε_ijk

where:
  y_ijk = observed phenotype
  g_i(jk) = area-specific genetic value (breeding value)
  e_j = environment main effect
  b_k(j) = area effect nested in environment
  ε_ijk = residual error
```

**Variance components**:
```
Var(y_ijk) = σ²_g(j) + σ²_env + σ²_area + σ²_ε
           = σ²_g(j) + 2σ²_g(j) + 0.5σ²_g(j) + σ²_ε
           = 3.5σ²_g(j) + σ²_ε
```

### Heritability Specification

**Treatment parameter**: h² = 0.2, 0.4, 0.6, 0.8

**Definition**: Narrow-sense heritability on a **per-plot basis**:

```
h²_j = σ²_g(j) / (σ²_g(j) + σ²_ε)
```

**Error variance derived from target heritability**:

```r
σ²_ε = ((1 - h²) / h²) × σ²_g(j)
```

**Example calculations**:

| h² | σ²_g | σ²_ε | Total phenotypic variance |
|----|------|------|---------------------------|
| 0.2 | 1.0 | 4.00 | 5.00 |
| 0.4 | 1.0 | 1.50 | 2.50 |
| 0.6 | 1.0 | 0.67 | 1.67 |
| 0.8 | 1.0 | 0.25 | 1.25 |

**Note**: This is plot-basis heritability, not entry-mean heritability. Entry-mean heritability would be higher for replicated entries.

### Residual Error Simulation

**Two-step process to ensure exact variance**:

```r
# Step 1: Sample from standard normal
ε_raw ~ N(0, 1)

# Step 2: Standardize and rescale
ε_standardized = (ε_raw - mean(ε_raw)) / sd(ε_raw)
ε_ijk = ε_standardized × sqrt(σ²_ε)
```

**Rationale**:
- Simple sampling from N(0, σ²_ε) produces stochastic variance
- Standardization ensures **exact** error variance in each simulation
- Eliminates variance in error variance across simulation replicates
- Critical for fair comparison of designs across replicates

**Verification**:
```r
var(ε_ijk)  # Always exactly equals σ²_ε
```

### Phenotype Generation Algorithm

**Step-by-step process for each simulation run**:

```r
for (each combination of h², r_{k_{1}k_{2}}, replication level) {
  
  # 1. Generate genetic covariance matrix
  Σ_g <- kronecker(J_macro, K_intra)
  
  # 2. Sample marker effects
  β <- mvrnorm(n_markers, mu = 0, Sigma = Σ_g)
  
  # 3. Calculate area-specific breeding values
  BV <- geno %*% β
  
  # 4. Calculate true breeding value (mean across areas)
  trueBV <- rowMeans(BV)
  
  # 5. For each environment:
  for (j in 1:n_environments) {
    
    # Calculate genetic variance in this environment
    σ²_g_j <- var(BV[, areas_in_env_j])
    
    # Sample environment effect
    e_j ~ N(0, 2 × σ²_g_j)
    
    # For each area in environment:
    for (k in 1:n_areas_per_env) {
      
      # Sample area effect
      b_k_j ~ N(0, 0.5 × σ²_g_j)
      
      # Calculate target error variance
      σ²_ε <- ((1 - h²) / h²) × σ²_g_j
      
      # Sample and standardize residual error
      ε_raw ~ N(0, 1)
      ε <- (ε_raw - mean(ε_raw)) / sd(ε_raw) × sqrt(σ²_ε)
      
      # Generate phenotypes
      y <- BV[, area_k_in_env_j] + e_j + b_k_j + ε
    }
  }
}
```

### Quality Control Checks

After phenotype generation, verify:

```r
# 1. Achieved heritability matches target
h²_realized <- var(BV) / var(y)
stopifnot(abs(h²_realized - h²_target) < 0.01)

# 2. Error variance matches target
σ²_ε_realized <- var(ε)
stopifnot(abs(σ²_ε_realized - σ²_ε_target) < 0.001)

# 3. Genetic correlations match target
cor_realized <- cor(BV[, area1], BV[, area2])
stopifnot(abs(cor_realized - r_{k_{1}k_{2}}_target) < 0.05)

# 4. No missing values generated
stopifnot(!any(is.na(y)))
```

---

## Statistical Models

### Model 1: Conventional Selection (BLUP)

**Purpose**: Estimate breeding values using phenotypic data only (no genomic information)

#### Model Specification

**Fixed effects**:
```
y_ijk = μ + e_i + ...
```

**Random effects**:
```
... + b_k(i) + g_j(i) + ε_ijk

where:
  b_k(i) ~ N(0, B)    # Area effects, environment-specific variances
  g_j(i) ~ N(0, G_e)  # Genotype effects, environment-specific variances
  ε_ijk ~ N(0, R)     # Residual errors, environment-specific variances
```

#### Variance Structure

**Block variance (B)**:
```
B = diag(σ²_b1, σ²_b2, ..., σ²_be)
```
- Diagonal matrix
- Environment-specific block variances estimated separately

**Genetic variance (G_e)**:
```
G_e = diag(σ²_g1, σ²_g2, ..., σ²_ge)
```
- Diagonal matrix
- Environment-specific genetic variances estimated separately
- Lines assumed **unrelated** (G_m = I)

**Residual variance (R)**:
```
R = diag(σ²_ε1, σ²_ε2, ..., σ²_εe) ⊗ I_n
```
- Diagonal for environments
- Identity for observations within environment

#### ASReml-R Implementation

**For high/intermediate replication** (≥2 reps):
```r
BLUP <- asreml(
  fixed = phenotype ~ location,
  random = ~ location:block + at(location):idv(germplasmName),
  residual = ~ dsum(~ units | location),
  data = breedSim,
  na.action = na.method(y = "omit", x = "omit")
)
```

**For low replication** (unreplicated entries):
```r
BLUP <- asreml(
  fixed = phenotype ~ location,
  random = ~ at(location):idv(germplasmName),  # No block term
  residual = ~ dsum(~ units | location),
  data = breedSim,
  na.action = na.method(y = "omit", x = "omit")
)
```

**Key model components**:
- `at(location):idv(germplasmName)`: Heterogeneous genetic variances by environment
- `location:block`: Block effects nested in environment
- `dsum(~ units | location)`: Heterogeneous residual variances by environment
- `idv()`: Identity variance structure (lines unrelated)

#### Model Fitting Strategy

**For RR designs**: Fit separate models for each trial (preliminary, advanced)
```r
# Preliminary trial model
BLUPprelim <- asreml(...)

# Advanced trial model
BLUPadv <- asreml(...)

# Combine predictions
BLUPpred <- bind_rows(BLUPprelim$predictions, BLUPadv$predictions)
```

**For CR designs**: Fit single model across all cohorts
```r
BLUP <- asreml(...)
BLUPpred <- BLUP$predictions
```

#### Prediction Accuracy

```r
# Extract predicted breeding values
BV_predicted <- predict(BLUP, classify = "germplasmName")$pvals$predicted.value

# Calculate accuracy
r_bv <- cor(BV_predicted, trueBV, method = "pearson")
```

**Reported metrics**:
- **Overall**: Correlation across all experimental lines
- **By cohort**: Separate correlations for S1, S2, S3, S4
- **By trial** (RR only): Separate correlations for preliminary and advanced

### Model 2: Genomic Selection (GBLUP)

**Purpose**: Estimate breeding values using both phenotypic and genomic relationship information

#### Model Specification

Identical structure to BLUP, but with genomic relationship matrix:

**Random effects**:
```
g_j(i) ~ N(0, G_e ⊗ G_m)

where:
  G_e = diag(σ²_g1, ..., σ²_ge)  # Environment-specific genetic variances
  G_m = A                         # Genomic relationship matrix
```

#### Genomic Relationship Matrix

**Calculation** (VanRaden's method):
```r
# Using rrBLUP package
library(rrBLUP)
K2 <- A.mat(geno - 1)  # Subtract 1 to center markers

# Manual calculation
M_centered <- geno - matrix(colMeans(geno), nrow(geno), ncol(geno), byrow = TRUE)
K2 <- (M_centered %*% t(M_centered)) / (2 × sum(p × (1 - p)))
```

**Properties**:
- Symmetric positive semi-definite matrix
- Diagonal elements ≈ 1 + inbreeding coefficient
- Off-diagonal elements = genomic relationships between lines
- Captures realized genomic similarity more precisely than pedigree

**Matrix dimensions**: n_lines × n_lines (4,102 × 4,102)

#### ASReml-R Implementation

**All replication levels**:
```r
GBLUP <- asreml(
  fixed = phenotype ~ location,
  random = ~ location:block + at(location):vm(germplasmName, K2),
  residual = ~ dsum(~ units | location),
  data = breedSim,
  na.action = na.method(y = "omit", x = "omit")
)
```

**Key differences from BLUP**:
- `vm(germplasmName, K2)`: Variance model using genomic relationship matrix
- Lines now **related** through realized genomic similarity
- Single model fits both RR and CR designs (no separate trial models needed)

#### Convergence Settings

GBLUP requires more computational resources:
```r
asreml.options(
  pworkspace = "8gb",      # Increased workspace
  ai.sing = TRUE,          # Allow singular AI matrix
  fail = "soft",           # Continue with warnings
  threads = 6              # Parallel processing
)
```

#### Prediction Accuracy

```r
# Extract genomic estimated breeding values (GEBVs)
GEBV <- predict(GBLUP, classify = "germplasmName")$pvals$predicted.value

# Calculate accuracy
r_bv <- cor(GEBV, trueBV, method = "pearson")
```

**Reported metrics**:
- **Overall**: Correlation across all experimental lines
- **By cohort**: Separate correlations for S1, S2, S3, S4

**Note**: No separate trial-level accuracies reported for GBLUP since genomic relationships connect all lines regardless of design.

### Model 3: Sparse Testing

**Purpose**: Evaluate genomic prediction when not all lines are tested in all environments

#### Sparse Testing Scenario

**Training set**:
- S3 and S4 cohorts tested in **all 5 environments**
- Provides phenotypic data for model training

**Test set**:
- S1 and S2 cohorts tested in **1 of 5 environments only**
- Requires prediction for 4 untested environments

**Two sparse testing designs evaluated**:

1. **Sparse (one location)**:
   ```
   S1, S2: Only environment 1 has phenotypes
   S3, S4: All 5 environments have phenotypes
   ```

2. **Sparse (no preliminary)**:
   ```
   S1, S2: No phenotypes in any environment
   S3, S4: All 5 environments have phenotypes
   ```

#### Data Structure

```r
# Sparse testing with one location
breedSim1Loc <- breedSim %>%
  mutate(GSset = if_else(test == "prelim" & location != 1, "test", "train")) %>%
  mutate(phenotype = if_else(GSset == "test", NA_real_, phenotype))

# Sparse testing with no preliminary phenotypes
breedSimNoLoc <- breedSim %>%
  mutate(GSset = if_else(test == "prelim", "test", "train")) %>%
  mutate(phenotype = if_else(GSset == "test", NA_real_, phenotype))
```

#### ASReml-R Implementation

**Critical**: Use `na.action = na.method(y = "include")` to retain lines with missing phenotypes

```r
Sparse <- asreml(
  fixed = phenotype ~ location,
  random = ~ location:block + at(location):vm(germplasmName, K2),
  residual = ~ dsum(~ units | location),
  data = breedSim1Loc,
  na.action = na.method(y = "include", x = "include")  # Include missing y
)
```

**Model behavior**:
- Lines with phenotypes contribute to variance component estimation
- Lines without phenotypes receive predictions based on:
  - Genomic relationships with phenotyped lines
  - Environment main effects
  - Estimated variance components

#### Prediction Accuracy

**Only evaluate test set**:
```r
# Extract predictions for test set only
resultsByEntry <- breedSim1Loc %>%
  filter(GSset == "test" & cohort != "check") %>%
  distinct(germplasmName, .keep_all = TRUE) %>%
  left_join(SparsePred, by = "germplasmName")

# Calculate accuracy for untested lines
r_bv_sparse <- cor(resultsByEntry$predicted.value, 
                   resultsByEntry$trueBV, 
                   method = "pearson")
```

**Reported metrics**:
- **Test set overall**: Correlation for S1 and S2 combined
- **By cohort**: Separate correlations for S1 and S2

**Interpretation**: These accuracies represent the program's ability to predict performance in untested environments using genomic relationships.

### Variance Component Estimation

All models estimate environment-specific variance components:

**Genetic variance by environment**:
```
σ²_g(j) for j = 1, ..., 5
```

**Block variance by environment**:
```
σ²_b(j) for j = 1, ..., 5
```

**Residual variance by environment**:
```
σ²_ε(j) for j = 1, ..., 5
```

**Heritability by environment** (calculated post-hoc):
```
h²_j = σ²_g(j) / (σ²_g(j) + σ²_ε(j) / n_rep(j))
```

---

## Simulation Scenarios

### Complete Factorial Design

**Total unique scenarios**:

| Model Type | h² levels | r_{k_{1}k_{2}} levels | Rep levels | Design types | Total |
|------------|-----------|----------------------|------------|--------------|-------|
| BLUP | 4 | 5 | 3 | 2 | 120 |
| GBLUP | 4 | 5 | 1 (low) | 2 | 40 |
| Sparse | 4 | 5 | 1 (low) | 2 | 40 |

### Scenario Specifications

#### BLUP Scenarios (120 total)

**Fixed parameters**:
- Environments: 5
- Macro-GxE: r = 0.5
- Design types: RR, CR

**Varied parameters**:
```
h² = {0.2, 0.4, 0.6, 0.8}                    (4 levels)
r_{k_{1}k_{2}} = {0.2, 0.4, 0.6, 0.8, 1.0}   (5 levels)
Replication = {Low, Intermediate, High}       (3 levels)
Design = {RR, CR}                             (2 levels)

Total: 4 × 5 × 3 × 2 = 120 scenarios
```

**Simulation replicates**: 500 per scenario

**Rationale for high replication**:
- BLUP computationally inexpensive (no genomic matrix)
- High replication reduces Monte Carlo error
- Enables precise estimation of small effect sizes

**Total BLUP model fits**: 120 scenarios × 500 reps = 60,000

#### GBLUP Scenarios (40 total)

**Fixed parameters**:
- Environments: 5
- Macro-GxE: r = 0.5
- Replication: Low only
- Design types: RR, CR

**Varied parameters**:
```
h² = {0.2, 0.4, 0.6, 0.8}                    (4 levels)
r_{k_{1}k_{2}} = {0.2, 0.4, 0.6, 0.8, 1.0}   (5 levels)
Design = {RR, CR}                             (2 levels)

Total: 4 × 5 × 2 = 40 scenarios
```

**Simulation replicates**: 50 per scenario

**Rationale for reduced scope**:
- GBLUP computationally intensive
- Initial BLUP analysis showed no replication × design interaction
- Focus on most relevant breeding scenarios (low replication with genomic data)

**Total GBLUP model fits**: 40 scenarios × 50 reps = 2,000

#### Sparse Testing Scenarios (40 total)

**Fixed parameters**:
- Environments: 5 (S1/S2 tested in 1 only)
- Macro-GxE: r = 0.5
- Replication: Low only
- Design types: RR, CR

**Varied parameters**:
```
h² = {0.2, 0.4, 0.6, 0.8}                    (4 levels)
r_{k_{1}k_{2}} = {0.2, 0.4, 0.6, 0.8, 1.0}   (5 levels)
Design = {RR, CR}                             (2 levels)
Model = {Sparse1Loc, SparseNoPrelim}          (2 models per scenario)

Total scenarios: 4 × 5 × 2 = 40
Total model fits: 40 × 2 models = 80
```

**Simulation replicates**: 100 per scenario

**Rationale for increased replication**:
- Sparse testing has higher Monte Carlo variance
- Smaller test set (S1, S2 only) requires more replicates for stable estimates
- Each replicate fits 2 models (Sparse1Loc vs SparseNoPrelim)

**Total sparse testing model fits**: 40 scenarios × 100 reps × 2 models = 8,000

### Computational Resource Requirements

#### Hardware Specifications

**Recommended minimum**:
```
CPU: 6+ cores
RAM: 32 GB
Storage: 100 GB for outputs
OS: Linux/Unix preferred (for parallel processing)
```

#### Execution Time Estimates

**BLUP simulations**:
```
Per scenario: ~45 minutes (500 replicates)
Total time: 120 scenarios × 45 min ≈ 90 hours ≈ 3.75 days
Parallelization: 6 cores reduces to ~15 hours
```

**GBLUP simulations**:
```
Per scenario: ~30 minutes (50 replicates)
Total time: 40 scenarios × 30 min ≈ 20 hours ≈ 1 day
Parallelization: 4 cores reduces to ~5 hours
```

**Sparse testing simulations**:
```
Per scenario: ~45 minutes (100 replicates × 2 models)
Total time: 40 scenarios × 45 min ≈ 30 hours ≈ 1.25 days
Parallelization: 4 cores reduces to ~7.5 hours
```

**Total project time**: ~5 days on single workstation, ~1.5 days with parallelization

#### Memory Usage

**Peak memory per model fit**:
```
BLUP: ~2 GB
GBLUP: ~8 GB
Sparse: ~8 GB
```

**Parallel processing considerations**:
- BLUP: Can run 6+ simultaneous jobs (2 GB × 6 = 12 GB)
- GBLUP: Limited to 4 simultaneous jobs (8 GB × 4 = 32 GB)
- Sparse: Limited to 4 simultaneous jobs (8 GB × 4 = 32 GB)

### Error Handling

**ASReml convergence failures**:
- Tracked separately using `safely()` wrapper
- Saved to error log files
- Typical failure rate: <1% of model fits
- Common causes:
  - Extremely low heritability (h² = 0.2) with low r_{k_{1}k_{2}} (0.2)
  - Singularity in random effects structure
  - Insufficient replication for variance component estimation

**Quality control**:
```r
# Track successful vs failed fits
allOutput <- simulations %>% future_map(gsAccSafe)
results <- allOutput %>% future_map("result") %>% compact()
errors <- allOutput %>% future_map("error")

# Save separately
saveRDS(errors, "errors_log.Rdata")
write_csv(results, "results.csv")
```

---

## Output Specifications

### File Organization

```
output/
├── BLUPoutput/
│   ├── BLUPresults/
│   │   ├── results_1_h0.2_1-16.csv
│   │   ├── results_1_h0.4_1-16.csv
│   │   └── ... (450 iterations × 4 h² levels)
│   ├── BLUPerrors/
│   │   ├── errors_1_h0.2_1-16.Rdata
│   │   └── ...
│   └── BLUPfinal/
│       └── breedSimResult_BLUP_1-16.csv
├── gBLUPoutput/
│   ├── gBLUPresults/
│   │   ├── gBLUPresultsV2_12-21_1.csv
│   │   └── ... (50 iterations)
│   ├── gBLUPerrors/
│   └── gBLUPresults_final/
│       └── breedSimResult_gBLUP_allIterationsV2_12-21.csv
└── SparseOutput/
    ├── SparseResults/
    │   ├── SparseResults_1_h0.2_1-16.csv
    │   └── ... (100 iterations × 4 h² levels)
    ├── SparseErrors/
    └── SparseResults_final/
        └── breedSimResult_Sparse_1-16.csv
```

### Output Data Structure

Each results CSV contains the following columns:

#### Core Identification Variables

| Column | Type | Description | Example Values |
|--------|------|-------------|----------------|
| `iteration` | integer | Simulation replicate number | 1, 2, ..., 500 |
| `design` | character | Randomization scheme | "RCBD", "PREP" |
| `heritability` | numeric | Target heritability | 0.2, 0.4, 0.6, 0.8 |
| `nLoc` | integer | Number of environments | 5 |
| `macroGxE` | numeric | Inter-environment correlation | 0.5 |
| `microGxE` | numeric | Intra-environment correlation (r_{k_{1}k_{2}}) | 0.2, 0.4, 0.6, 0.8, 1.0 |
| `repCat` | character | Replication level | "1-1-2-2", "2-2-2-2", "2-2-2-3" |

#### Accuracy Metrics

| Column | Type | Description | Values |
|--------|------|-------------|--------|
| `cor` | numeric | Prediction accuracy (Pearson r) | -1 to 1 |
| `group` | character | Subset for accuracy calculation | "overall", "S1", "S2", "S3", "S4", "prelim", "adv", "testSet" |
| `model` | character | Prediction method | "BLUP", "gBLUP", "Sparse", "SparseNoPrelim" |

#### Optional Metadata

| Column | Type | Description |
|--------|------|-------------|
| `time` | POSIXct | Model completion timestamp |

### Replication Category Encoding

The `repCat` column encodes replication structure:

| Code | S1 reps | S2 reps | S3 reps | S4 reps | Description |
|------|---------|---------|---------|---------|-------------|
| `1-1-2-2` | 1 | 1 | 2 | 2 | Low replication |
| `2-2-2-2` | 2 | 2 | 2 | 2 | Intermediate replication |
| `2-2-2-3` | 2 | 2 | 2 | 3 | High replication |

### Accuracy Group Definitions

| Group Value | Description | Used In |
|-------------|-------------|---------|
| `overall` | All experimental lines (S1-S4 combined) | BLUP, GBLUP |
| `S1` | Stage 1 cohort only | All models |
| `S2` | Stage 2 cohort only | All models |
| `S3` | Stage 3 cohort only | BLUP, GBLUP |
| `S4` | Stage 4 cohort only | BLUP, GBLUP |
| `prelim` | Preliminary trial (S1 + S2) | BLUP (RR only) |
| `adv` | Advanced trial (S3 + S4) | BLUP (RR only) |
| `testSet` | Untested lines (S1 + S2 in sparse testing) | Sparse only |

### Example Output Records

**BLUP result (RR design)**:
```csv
iteration,design,heritability,nLoc,macroGxE,microGxE,repCat,group,cor,model
1,RCBD,0.6,5,0.5,0.8,2-2-2-2,overall,0.847,BLUP
1,RCBD,0.6,5,0.5,0.8,2-2-2-2,S1,0.823,BLUP
1,RCBD,0.6,5,0.5,0.8,2-2-2-2,prelim,0.831,BLUP
```

**GBLUP result (CR design)**:
```csv
iteration,design,heritability,nLoc,macroGxE,microGxE,repCat,group,cor,model,time
1,PREP,0.4,5,0.5,0.6,1-1-2-2,overall,0.888,gBLUP,2025-01-15 14:32:18
1,PREP,0.4,5,0.5,0.6,1-1-2-2,S1,0.862,gBLUP,2025-01-15 14:32:18
```

**Sparse testing result**:
```csv
iteration,design,heritability,nLoc,macroGxE,microGxE,repCat,group,cor,model
1,PREP,0.8,5,0.5,1.0,1-1-2-2,testSet,0.812,Sparse
1,PREP,0.8,5,0.5,1.0,1-1-2-2,S1,0.798,Sparse
1,PREP,0.8,5,0.5,1.0,1-1-2-2,testSet,0.754,SparseNoPrelim
```

### Data Aggregation

For analysis, data typically aggregated by averaging across iterations:

```r
# Calculate mean accuracy and standard error by scenario
results_summary <- results %>%
  group_by(design, heritability, microGxE, repCat, group, model) %>%
  summarize(
    mean_cor = mean(cor, na.rm = TRUE),
    se_cor = sd(cor, na.rm = TRUE) / sqrt(n()),
    n_iterations = n()
  )
```

---

## Parameter Sensitivity

### Key Findings from Sensitivity Analysis

#### 1. Heritability (h²)

**Effect magnitude**: Strongest main effect on prediction accuracy

**Pattern**:
```
r_bv increases monotonically with h²:
  h² = 0.2: r_bv ≈ 0.55-0.65
  h² = 0.4: r_bv ≈ 0.70-0.80
  h² = 0.6: r_bv ≈ 0.80-0.88
  h² = 0.8: r_bv ≈ 0.88-0.95
```

**Design interaction**: **Non-significant**
- RR and CR respond nearly identically to changes in h²
- Difference-in-differences (DiD) coefficient: δ̂ = 0.001, p > 0.05
- Both designs suffer equally under low heritability conditions

**Biological interpretation**:
- Low h² means phenotype dominated by non-genetic variation
- Both designs struggle when signal-to-noise ratio is poor
- Design choice doesn't compensate for fundamental measurement limitations

**Recommendation**: Heritability can be treated as **nuisance parameter** when comparing designs

#### 2. Intra-Environment Genetic Correlation (r_{k_{1}k_{2}})

**Effect magnitude**: **Critical driver of design performance differences**

**Pattern**:
```
Design advantage increases as r_{k_{1}k_{2}} decreases:
  r_{k_{1}k_{2}} = 1.0: RR ≈ CR (identical performance)
  r_{k_{1}k_{2}} = 0.8: CR advantage = +2.5%
  r_{k_{1}k_{2}} = 0.6: CR advantage = +8.8%
  r_{k_{1}k_{2}} = 0.4: CR advantage = +18.0%
  r_{k_{1}k_{2}} = 0.2: CR advantage = +18.6%
```

**Design interaction**: **Highly significant**
- BLUP DiD: δ̂ = 0.082, p < 0.001
- GBLUP DiD: δ̂ = 0.003, p > 0.10 (non-significant)
- Sparse DiD: δ̂ = 0.018, p < 0.01

**Biological interpretation**:
- r_{k_{1}k_{2}} represents consistency of genetic rankings across spatial units (trials or blocks)
- Low r_{k_{1}k_{2}} = strong localized genotype-by-environment interaction
- When ranks change substantially across space, design matters more

**Mechanism**:
```
In RR designs (low r_{k_{1}k_{2}}):
  - Cohorts spatially segregated
  - Genetic differences confounded with trial effects
  - Variance partitioning problematic
  → Inflated or deflated breeding value estimates

In CR designs (low r_{k_{1}k_{2}}):
  - Cohorts interspersed within same space
  - Trial effects orthogonal to genetic effects
  - Clean variance partitioning
  → Unbiased breeding value estimates
```

**Recommendation**: r_{k_{1}k_{2}} is the **primary consideration** for design choice

#### 3. Replication Level

**Effect magnitude**: Moderate main effect; small design interaction

**Pattern**:
```
Accuracy improves with replication (logarithmic relationship):
  Low → Intermediate: Large gain (+8-10%)
  Intermediate → High: Modest gain (+2-3%)
```

**Design interaction**: **Significant but small**
- DiD coefficient: δ̂ = 0.005, p < 0.001
- Differential effect emerges only at low replication

**Design-specific patterns**:

| Replication | RR accuracy loss | CR accuracy loss |
|-------------|------------------|------------------|
| High → Intermediate | -8.2% | -8.0% |
| Intermediate → Low | -2.98% | -1.94% |

**Biological interpretation**:
- First replicate provides largest information gain (Yan 2021)
- Unreplicated entries require spatial or genetic borrowing
- CR better leverages replicated entries to aid unreplicated entries

**Mechanism (low replication)**:
```
In RR with unreplicated entries:
  - Requires checks within each trial
  - Estimates trial-specific error variance
  - Spatial borrowing limited to within-trial

In CR with unreplicated entries (p-rep):
  - Replicated S3/S4 serve as "checks"
  - Unified error variance across all cohorts
  - Spatial borrowing across entire trial
```

**Recommendation**: Design choice most critical when **low replication unavoidable**

---

## Appendices

### Appendix A: Statistical Notation Guide

| Symbol | Description |
|--------|-------------|
| **Indices** | |
| i | Line/genotype index (1, ..., n) |
| j | Environment index (1, ..., e) |
| k | Area/block index within environment (1, ..., a) |
| l | Marker index (1, ..., m) |
| **Observed Variables** | |
| y_ijk | Phenotype of line i in environment j, area k |
| M | Marker matrix (n × m) |
| **Genetic Parameters** | |
| g_i(jk) | Breeding value of line i in environment j, area k |
| trueBV_i | True breeding value (mean across areas) |
| β_l | Effect of marker l |
| σ²_g | Genetic variance |
| σ²_GxE | Genotype-by-environment variance |
| r_{k_{1}k_{2}} | Genetic correlation between areas (intra-environment) |
| r | Genetic correlation between environments (inter-environment) |
| h² | Narrow-sense heritability (plot basis) |
| **Environmental Parameters** | |
| e_j | Main effect of environment j |
| b_k(j) | Main effect of area k within environment j |
| ε_ijk | Residual error |
| σ²_env | Environment variance |
| σ²_area | Area/block variance |
| σ²_ε | Residual error variance |
| **Matrix Structures** | |
| J | Macro-GxE covariance matrix (e × e) |
| K | Intra-environment GxE covariance matrix (a × a) |
| Σ_g | Full genetic covariance matrix = J ⊗ K |
| K2 (or A) | Genomic relationship matrix (n × n) |
| G_e | Environment-specific genetic variance matrix |
| G_m | Marker-based relationship matrix |
| B | Block variance matrix |
| R | Residual variance matrix |
| **Accuracy Metrics** | |
| r_bv | Prediction accuracy (correlation) |
| r̄_bv | Mean prediction accuracy across iterations |
| **DiD Parameters** | |
| Z | Randomization scheme (0=RR, 1=CR) |
| T | Parameter condition (0=optimal, 1=suboptimal) |
| δ̂_DiD | Difference-in-differences coefficient |

### Appendix B: ASReml-R Function Reference

#### Key Functions Used

**Model fitting**:
```r
asreml(fixed, random, residual, data, na.action, ...)
```

**Variance structures**:
```r
at(factor):term         # Heterogeneous variance by factor levels
idv(germplasmName)      # Identity variance (unrelated)
vm(germplasmName, K2)   # Variance model with relationship matrix
dsum(~ units | factor)  # Diagonal (heterogeneous) residual by factor
```

**Prediction**:
```r
predict.asreml(object, classify, ...)
```

**Options**:
```r
asreml.options(
  pworkspace = "8gb",   # Workspace memory
  ai.sing = TRUE,       # Allow singular AI matrix
  fail = "soft",        # Continue with warnings
  threads = 6           # Parallel threads
)
```

### Appendix C: Software Versions

**Core dependencies**:
```
R version: 4.2.1 (2022-06-23)
ASReml-R: 4.2.0.290 (commercial license required)
Platform: x86_64-pc-linux-gnu (64-bit)
```

**Key packages**:
```
data.table_1.14.8
tidyverse_2.0.0
  ├─ dplyr_1.1.2
  ├─ tidyr_1.3.0
  ├─ ggplot2_3.4.2
  └─ purrr_1.0.1
magrittr_2.0.3
furrr_0.3.1
future_1.32.0
MASS_7.3-60
matrixcalc_1.0-6
MBESS_4.9.2
caret_6.0-94
```

**Session info for reproducibility**:
```r
sessionInfo()
# R version 4.2.1 (2022-06-23)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 22.04.2 LTS
```

### Appendix D: Troubleshooting Guide

#### Common Errors and Solutions

**1. ASReml License Issues**

```
Error: ASReml license not found or expired
```

**Solution**:
```r
# Check license status
asreml::asreml.license.status()

# Set license file location (if needed)
options(asreml.license.file = "/path/to/asreml.lic")
```

**2. Memory Errors**

```
Error: cannot allocate vector of size X Gb
```

**Solution**:
```r
# Increase R memory limit (Windows)
memory.limit(size = 32000)

# Increase workspace for ASReml
asreml.options(pworkspace = "16gb")

# Reduce parallel workers
plan(multicore, workers = 2)  # Instead of 6
```

**3. Matrix Singularity**

```
Warning: Singular information matrix
```

**Causes**:
- Insufficient replication for variance component estimation
- Too many variance components relative to data
- Linear dependencies in fixed effects

**Solutions**:
```r
# Allow singular matrices
asreml.options(ai.sing = TRUE)

# Simplify model
# Remove block effects if unreplicated
model <- asreml(
  fixed = phenotype ~ location,
  random = ~ at(location):idv(germplasmName),  # Block term removed
  ...
)

# Check for sufficient data
table(data$location, data$germplasmName)  # Should have >1 obs per cell
```

**4. Convergence Failures**

```
Warning: Model did not converge
```

**Solutions**:
```r
# Increase iteration limit
asreml.options(maxiter = 50)  # Default is 20

# Use different starting values
model <- update(model, start.values = TRUE)

# Simplify variance structure
# Use homogeneous residual variance instead of heterogeneous
residual = ~ units  # Instead of dsum(~ units | location)
```

**5. Genomic Relationship Matrix Issues**

```
Error: K2 matrix not positive definite
```

**Solutions**:
```r
# Check matrix properties
isSymmetric(K2)
min(eigen(K2)$values)  # Should be > 0

# Correct if needed
library(Matrix)
K2_corrected <- as.matrix(nearPD(K2, keepDiag = TRUE)$mat)

# Add small value to diagonal (regularization)
K2_reg <- K2 + diag(0.001, nrow(K2))
```

**6. Factor Level Mismatches**

```
Error: germplasmName levels don't match between data and K2
```

**Solutions**:
```r
# Ensure exact matching
stopifnot(all(rownames(K2) == levels(data$germplasmName)))

# Reorder if needed
K2 <- K2[levels(data$germplasmName), levels(data$germplasmName)]

# Or use factor levels from data
data$germplasmName <- factor(data$germplasmName, 
                             levels = rownames(K2))
```

**7. Parallel Processing Errors**

```
Error: node produced an error
```

**Solutions**:
```r
# Use safer error handling
gsAccSafe <- safely(gsAcc)
results <- map(scenarios, gsAccSafe)

# Check for errors
errors <- results %>% map("error") %>% compact()
print(errors)

# Reduce workers if memory constrained
plan(multicore, workers = 2)

# Or use sequential processing for debugging
plan(sequential)
```

#### Data Validation Checks

**Pre-simulation checks**:
```r
# 1. Marker data quality
stopifnot(!any(is.na(geno)))
stopifnot(all(geno %in% c(-1, 0, 1)))
print(paste("MAF range:", 
            round(min(colMeans(geno + 1) / 2), 3), "to",
            round(max(colMeans(geno + 1) / 2), 3)))

# 2. Genomic relationship matrix
stopifnot(isSymmetric(K2))
stopifnot(min(eigen(K2)$values) > -1e-8)
print(paste("Relationship range:", 
            round(min(K2[lower.tri(K2)]), 3), "to",
            round(max(K2[lower.tri(K2)]), 3)))

# 3. Cohort structure
cohort_sizes <- table(uniqueLines$cohort)
print(cohort_sizes)
stopifnot(all(c("S1", "S2", "S3", "S4") %in% names(cohort_sizes)))
```

**Post-simulation checks**:
```r
# 1. Achieved heritability
h2_check <- results %>%
  group_by(heritability) %>%
  summarize(mean_accuracy = mean(cor))
print(h2_check)
# Should show monotonic increase

# 2. Correlation range
cor_range <- range(results$cor, na.rm = TRUE)
stopifnot(cor_range[1] >= 0 & cor_range[2] <= 1)

# 3. Missing values
n_missing <- sum(is.na(results$cor))
prop_missing <- n_missing / nrow(results)
print(paste("Missing accuracies:", n_missing, 
            sprintf("(%.2f%%)", prop_missing * 100)))
stopifnot(prop_missing < 0.05)  # Less than 5% failure rate

# 4. Baseline comparison (r_{k_{1}k_{2}} = 1.0)
baseline <- results %>%
  filter(microGxE == 1.0) %>%
  group_by(design) %>%
  summarize(mean_cor = mean(cor))
print(baseline)
# RR and CR should be nearly identical at r_{k_{1}k_{2}} = 1.0
stopifnot(abs(diff(baseline$mean_cor)) < 0.01)
```

### Appendix E: Parameter Decision Tree

**Use this flowchart to select appropriate parameters for new simulations:**

```
START
│
├─ Q1: Do you have genotypic data?
│  ├─ YES: Proceed to Q2
│  └─ NO: Use BLUP only
│     └─ Include h² = {0.2, 0.4, 0.6, 0.8}
│
├─ Q2: What is your testing design?
│  ├─ Balanced MET (all lines in all environments)
│  │  └─ Use BLUP + GBLUP
│  │     ├─ Test replication levels: {Low, Int, High}
│  │     └─ Test r_{k_{1}k_{2}} = {0.2, 0.4, 0.6, 0.8, 1.0}
│  │
│  └─ Sparse testing (lines not in all environments)
│     └─ Use BLUP + GBLUP + Sparse
│        ├─ Focus on low replication
│        └─ Test r_{k_{1}k_{2}} = {0.2, 0.4, 0.6, 0.8, 1.0}
│
├─ Q3: How many environments realistic for your program?
│  ├─ 3-5: Use n = 5 (our default)
│  ├─ 6-10: Consider n = 7
│  └─ >10: Consider n = 10, reduce other parameters
│
├─ Q4: What is typical G×E in your system?
│  ├─ Strong (different regions/years): r = 0.3-0.5
│  ├─ Moderate (similar regions): r = 0.5-0.7
│  └─ Weak (controlled conditions): r = 0.7-0.9
│
└─ Q5: How many simulation replicates?
   ├─ Quick assessment: 50 iterations
   ├─ Publication quality: 100-500 iterations
   └─ High precision: 500-1000 iterations
      (depends on computational resources)
```

### Appendix F: Glossary of Terms

**Agricultural/Breeding Terms**:

| Term | Definition |
|------|------------|
| Breeding value | Genetic merit of individual for trait; sum of additive genetic effects |
| Cohort | Group of lines at same stage of testing/advancement |
| Entry | Line or variety included in yield trial |
| G×E interaction | Differential performance of genotypes across environments |
| Heritability (h²) | Proportion of phenotypic variance due to additive genetic effects |
| Multi-environment trial (MET) | Trials conducted across multiple location-year combinations |
| p-rep design | Partially replicated design with some entries unreplicated |
| Sparse testing | Testing strategy where not all lines evaluated in all environments |
| Stage of testing | Year/phase of advancement in breeding pipeline (S1, S2, S3, S4) |
| Yield trial | Field experiment to evaluate agronomic performance |

**Statistical Terms**:

| Term | Definition |
|------|------------|
| BLUP | Best Linear Unbiased Prediction |
| Complete randomization (CR) | All cohorts randomized together in single trial |
| DiD (Difference-in-Differences) | Causal inference method comparing treatment effects across groups |
| GBLUP | Genomic BLUP (uses genomic relationship matrix) |
| GRM | Genomic Relationship Matrix |
| IBD | Incomplete Block Design |
| RCBD | Randomized Complete Block Design |
| Restricted randomization (RR) | Cohorts randomized separately in different trials |
| r_bv | Prediction accuracy (correlation between predicted and true breeding values) |
| r_{k_{1}k_{2}} | Genetic correlation between areas within environments |

**Variance Components**:

| Symbol | Definition |
|--------|------------|
| σ²_g | Genetic variance |
| σ²_env | Environment main effect variance |
| σ²_area | Area/block effect variance |
| σ²_ε | Residual error variance |
| σ²_GxE | Genotype-by-environment interaction variance |

### Appendix G: Key Differences from Version 1.0

**Version 1.1 Updates**:

1. **Terminology Corrections**:
   - Changed "microGxE" → "r_{k_{1}k_{2}}" throughout document
   - Updated to match manuscript terminology: "genetic correlation between areas within environments"
   - Corrected variable names in code examples

2. **Model Name Changes**:
   - Changed "GEBV" → "Sparse Testing" or "Sparse" in all contexts
   - Updated "GEBVNoPrelim" → "SparseNoPrelim"
   - Revised section headers and table entries accordingly

3. **Output Structure Updates**:
   - Updated file paths: GEBVoutput → SparseOutput
   - Revised column names in example outputs
   - Updated accuracy group definitions

4. **Consistency Improvements**:
   - Aligned all parameter descriptions with manuscript V4
   - Standardized mathematical notation
   - Corrected cross-references throughout document

---

## References

**Key Methodological Citations**:

1. **Experimental Design**:
   - Clarke, G.P.Y. & Stefanova, K. (2011). Optimal design for early-generation plant-breeding trials with unreplicated or partially replicated test lines. *Australian & New Zealand Journal of Statistics*, 53(4), 461-480.
   - Cullis, B.R., Smith, A.B., & Coombes, N.E. (2006). On the design of early generation variety trials with correlated data. *Journal of Agricultural, Biological, and Environmental Statistics*, 11(4), 381-393.
   - Piepho, H.P. & Williams, E.R. (2006). A comparison of experimental designs for selection in breeding trials with nested treatment structure. *Theoretical and Applied Genetics*, 113, 1505-1513.

2. **Genomic Selection**:
   - Combs, E. & Bernardo, R. (2013). Accuracy of genomewide selection for different traits with constant population size, heritability, and number of markers. *The Plant Genome*, 6(1), 1-7.
   - Atanda, S.A., et al. (2022). Sparse testing using genomic prediction improves selection for breeding targets in elite spring wheat. *Theoretical and Applied Genetics*, 135, 1939-1950.
   - Endelman, J.B. (2011). Ridge regression and other kernels for genomic selection with R package rrBLUP. *The Plant Genome*, 4(3), 250-255.

3. **Statistical Methods**:
   - Smith, A.B., et al. (2007). Varietal selection for perennial crops where data relate to multiple harvests from a series of field trials. *Euphytica*, 157, 253-266.
   - Rothbard, S., et al. (2024). A tutorial on applying the difference-in-differences method to health data. *Current Epidemiology Reports*, 11(2), 85-95.
   - Gilmour, A.R., et al. (2018). ASReml User Guide Release 4.2. VSN International Ltd.

4. **Empirical Wheat Data**:
   - Lado, B., et al. (2016). Strategies for selecting crosses using genomic prediction in two wheat breeding programs. *The Plant Genome*, 9(2), 1-12.
   - Yan, W. (2021). Estimation of the optimal number of replicates in crop variety trials. *Frontiers in Plant Science*, 11, 590762.
   - Gaire, R., et al. (2021, 2022). Multi-trait genomic selection and fusarium head blight resistance in wheat. *The Plant Genome*.

---

## Document Information

**Document Title**: Simulation Parameter Specifications Guide

**Version**: 1.1

**Last Updated**: January 2025

**Authors**: 
- Arlyn J. Ackerman (Breeding Insight, Cornell University)
- Jessica Rutkoski (Department of Crop Sciences, University of Illinois at Urbana-Champaign)

**Corresponding Authors**:
- **Arlyn Ackerman**: aja258@cornell.edu
- **Jessica Rutkoski**: rutkoski@illinois.edu

**Version History**:
| Version | Date | Changes | Author |
|---------|------|---------|--------|
| 1.0 | January 2025 | Initial comprehensive guide | AJ Ackerman |
| 1.1 | January 2025 | Corrected r_{k_{1}k_{2}} terminology; changed GEBV to Sparse Testing | AJ Ackerman |

**Suggested Citation**:
Ackerman, A.J. & Rutkoski, J. (2025). Simulation Parameter Specifications Guide (Version 1.1). GitHub repository: [repository URL]

**License**: MIT License

**Repository**: [Insert GitHub repository URL]

---

**End of Simulation Parameter Specifications Guide v1.1**

This document provides comprehensive technical specifications for reproducing the simulation study described in "Randomization across breeding cohorts improves accuracy of conventional and genomic selection." For questions, clarifications, or to report errors, please contact the authors or open an issue on the GitHub repository.
