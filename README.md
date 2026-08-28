# Linear vs. Threshold Model GEBV Comparison for Categorical Traits

R scripts to compare (G)EBVs from linear and threshold models for binary and ordinal (polychotomous) traits, and to convert variance components between the liability and observed scales. Developed for the genomic evaluation of morphological traits in Nellore cattle using solutions from the BLUPF90 family of programs.

## Background

Threshold models are the standard for categorical traits but are computationally demanding. Linear models are a robust, faster alternative, yet their solutions are on the observed scale and must be approximated to the liability scale before conversion to probabilities.

For binary traits, Hidalgo et al. (2024a, 2024b) proposed and tested an approximation (AP1), and de Oliveira Padilha et al. (2026) derived an improved one (AP2) in which the observed-scale GEBV is rescaled by z, the standard normal density at the threshold:

GEBV_liability ≈ GEBV_observed / z

These scripts generalize AP2 to traits with more than two categories. The scaling factor z is the binary special case of the slope of the observed score on the standardized liability (Gianola, 1979):

b = Σₖ (cₖ − μY) · [φ(τₖ₋₁) − φ(τₖ)]

where cₖ are the category codes, τₖ the thresholds obtained from the cumulative category proportions, and φ the standard normal density. For K = 2, b reduces to z, so the same code handles binary and ordinal traits. Probabilities per category follow Gianola and Foulley (1983): Pₖ = Φ(τₖ − u) − Φ(τₖ₋₁ − u).

No variance components are required for the GEBV approximation: the liability scale is standardized (σ²ₑ = 1), and the thresholds carry all the population-level information needed. The approximation is a first-order linearization around the population mean and applies to any random effect in the model (e.g., permanent environment, contemporary group), since each observed-scale solution is the liability-scale effect multiplied by the same constant b.

## Scripts

### `gebv_comparison_polychotomous.R`

Compares GEBVs from linear and threshold models fitted to the same data.

Per trait, it:
1. Reads animal solutions (`effect == 2`) from both models and maps renumbered IDs back to original IDs via the `renadd*.ped` files
2. Approximates linear-model GEBVs to the liability scale (AP2 generalized, GEBV_o / b)
3. Converts both models' GEBVs to the probability of the best category and to the expected score E[Y|u]
4. Reports Spearman rank correlations on the observed, liability, probability, and expected-score scales, plus regression parameters
5. Saves scatter (probability and liability scales) and density plots as high-resolution TIFFs with transparent background

Usage: edit only the bottom section — the base path and one `run_trait()` call per trait with its directory, solution file names, and the vector of category counts from the data used in that analysis. Pass `codes =` if the category coding is not 1..K.

### `ordinal_convert.R`

Converts variance components from the liability to the observed scale for ordinal traits using the polychotomous slope method (Gianola, 1979), with an extension for repeatability models (permanent environmental variance) following the logic of Rutledge (1977). Each liability variance ratio maps to the observed scale as b² · (σ²ᵢ / σ²_total). Returns thresholds, observed-scale variance components, heritability, repeatability, and the residual breakdown into liability-mapped variance and discretization noise.

## Requirements

- R (≥ 4.0) with `dplyr`, `tidyr`, `ggplot2`, `readr`
- Solution files and `renadd*.ped` from BLUPF90+ / CBLUPF90 runs of the same model on the same data for both approaches

## Caveats

- The category counts supplied to each analysis must come from the data actually used in that model run, since thresholds and b depend on the category proportions.
- The scaling factor is fixed across the population; strong heterogeneity of category proportions across contemporary groups degrades the approximation (de Oliveira Padilha et al., 2026).
- For traits with extreme category concentration (analogous to extreme prevalence in binary traits), agreement between linear and threshold models is expected to weaken (Hidalgo et al., 2024b; de Oliveira Padilha et al., 2026).

## References

- de Oliveira Padilha, D. A., et al. 2026. Comparison of approximation methods for genomic estimated breeding values from observed to liability scales in dairy cattle health traits. J. Dairy Sci. 109:2787–2799. https://doi.org/10.3168/jds.2025-27502
- Dempster, E. R., and I. M. Lerner. 1950. Heritability of threshold characters. Genetics 35:212–236.
- Gianola, D. 1979. Heritability of polychotomous characters. Genetics 93:1051–1055.
- Gianola, D., and J. L. Foulley. 1983. Sire evaluation for ordered categorical data with a threshold model. Genet. Sel. Evol. 15:201–224.
- Hidalgo, J., et al. 2024a. Transforming estimated breeding values from observed to probability scale: how to make categorical data analyses more efficient. J. Anim. Sci. 102:skae307. https://doi.org/10.1093/jas/skae307
- Hidalgo, J., et al. 2024b. Converting estimated breeding values from the observed to probability scale for health traits. J. Dairy Sci. 107:9628–9637. https://doi.org/10.3168/jds.2024-24767
- Rutledge, J. J. 1977. Repeatability of threshold traits. Biometrics 33:395–399.
