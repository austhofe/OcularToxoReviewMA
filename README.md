DOI here when published

# About The Project

This is the research compendium for "Ocular Toxoplasmosis Infection Leading to Uveitis and or Lesions - A Systematic Review and Meta-Analysis" currently under submission. In this article we found 63 articles which describe uveitis and lesions following toxoplasmosis. 

PROSPERO Registration #: CRD42021289925

## Citation
Boyd K, Condrey K, Rosa Hernandez A, Austhof E, Lin L, Dehnbostel1 J, Hoffmann S, Flaxel C, Pogreba-Brown K. Ocular Toxoplasmosis Infection Leading to Uveitis and or Lesions - A Systematic Review and Meta-Analysis
Full Article link here when published

## Contents

This research compendium includes:
* Data extracted from included research articles
* Excel sheet with tables and summary measures
* Files from literature searches for each database
* Supplemental material file with additional tables and information

## Codes
The codes branch contains the codes for a Bayesian hierarchical model with a zero-inflated mixture prior, used in this project.

### Files

#### 1. `ma.prop.mix.R`

This script defines the primary modeling function: `ma.prop.mix()`, which performs a Bayesian meta-analysis of proportions using a **zero-inflated binomial mixture model** implemented in **JAGS**. The function is designed to handle datasets with potential structural zeros.

##### Function: `ma.prop.mix(e, n, data, mass0 = TRUE, ...)`

**Arguments**:
- `e`: Vector of event counts (number of cases).
- `n`: Vector of sample sizes.
- `data`: Optional data frame containing `e` and `n`.
- `mass0`: Logical; if `TRUE`, fits a zero-inflated mixture model assuming some studies contribute no signal ("mass at zero").
- Additional MCMC control arguments:
  - `n.adapt`, `n.chains`, `n.burnin`, `n.iter`, `n.thin`: Control JAGS sampling behavior.
  - `seed`: Sets the random seed for reproducibility.

**Workflow**:
1. Validates and cleans the input data, removing studies with missing entries.
2. Defines the JAGS model as a string, supporting both structural-zero and standard random-effects models.
3. Initializes JAGS chains with fixed random seeds.
4. Runs the JAGS model and extracts posterior draws for:
   - `p.mass`: Posterior mean probability of structural zero.
   - `prop`: Median prevalence estimate with a 95% credible interval.
   - `pred.nonmass` & `pred.overall`: Posterior predictive intervals.
5. Returns both raw posterior samples and summary statistics.

---

#### 2. `processing_code.Rmd`

This R Markdown document provides the full analysis pipeline for prevalent cases using the above model.

**Contents**:
- **Data Import**: Loads study-level data from an Excel file.
- **Data Cleaning**: Removes missing values and prepares inputs for modeling.
- **Model Fitting**: Applies `ma.prop.mix()` under both "mass at zero" and "non-mass" assumptions.
- **Posterior Analysis**: Extracts, summarizes, and interprets results.
- **Visualization** *(if applicable)*: Can include posterior distributions, predictive intervals, and summary tables.


### How to Run

1. Open R or RStudio and install any required packages:
   - `rjags`, `readxl`, `rmarkdown`, `ggplot2`, etc.

2. Source the model function:

```r
source("ma.prop.mix_updated seed.R")

## License

Distributed under the MIT License. See `LICENSE.txt` for more information.
Researchers interested in using this data for subsequent publications or research projects should contact the corresponding author, Kristen Pogreba-Brown, kpogreba@arizona.edu.


## Contact

Erika Austhof - barrette@arizona.edu

Project Link: ([OcularToxoReviewMA](https://github.com/austhofe/OcularToxoReviewMA/))
