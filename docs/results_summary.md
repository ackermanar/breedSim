# Research Results Summary

## Key Findings

This simulation study examined how experimental design affects prediction accuracy in wheat breeding programs, comparing **Complete Randomization (CR)** across breeding cohorts versus **Restricted Randomization (RR)** within separate trials.

### Main Results

**1. Design Performance Depends on Data Completeness**
- With limited phenotypic/genomic data, p-rep designs achieved **8.3 percentage points higher** prediction accuracy than RCBD designs (Δ% = +11.7%)
- Under low replication, the advantage increased to **10.1 percentage points** (Δ% = +15.7%)
- With complete data in GBLUP analysis, both designs performed similarly (r_bv = 0.888 vs 0.885)

**2. Genotype-by-Environment Interaction Drives Design Differences**
- Heritability showed strongest overall effect on prediction ability across all models
- However, GxE interaction produced the most significant design differences
- DiD analysis revealed highly significant interactions for BLUP (δ̂ᴰⁱᴰ = 0.078, p < .001) and sparse-teseted GBLUP (δ̂ᴰⁱᴰ = 0.018, p < .01)
- Restricted randomization can inflate genotype-environment covariance, creating prediction bias

**3. Practical Recommendations for Breeders**

Programs with comprehensive genomic/phenotypic data can use either design approach effectively. For resource-constrained programs or sparse testing scenarios, fully randomized designs (p-rep) provide meaningful advantages by:
- Better estimating environmental covariance
- Leveraging genomic relationships across environments
- Improving prediction stability under complex GxE interactions

For traits with low heritability, design choice becomes less critical as both approaches respond similarly.

---

**Simulation Parameters**: 5 environments, 4 breeding cohorts (S1-S4), heritability levels (0.2-0.8), genetic correlations (0.2-1.0), 500 iterations for BLUP, 50 iterations for genomic models.

**Analysis Methods**: Mixed models with environment-specific variance structures, difference-in-differences estimation for causal inference, multi-factor ANOVA for treatment effects.
