# A Bayesian Contamination Model for Serial Dilution Assays

Code for reproducing the results in the manuscript accepted for publication in *Biometrics*.

## Requirements

- R (version 4.0 or later)
- [RStan](https://mc-stan.org/rstan/) with a working C++ toolchain
- R packages: `tidyverse`, `rstan`, `ggpubr`, `psych`, `loo`, `shinystan`, `cowplot`, `Rmpfr`, `ggrepel`

## Repository structure

```
code/
  data_and_results/
    Generating results for figures.R   # fits the proposed model, saves model_fit.rds
    Bayesian_contamination_model.stan  # proposed model (Stan)
    derf.csv                           # standard + unknown assay data used in the analysis
    Documentation.txt
  figures/
    Figure 1.R ... Figure 5.R          # scripts that reproduce each figure
    Figure 4 Data.csv
    initial_model.stan                 # baseline model
    intermediate_model.stan
    readme.txt
  simulation/
    simulation_1/simulation_1_data_generation.R
    simulation_2/simulation_2_data_generation.R
data/
  Serial_Dilution_Assay_Sample_Data.csv
```

## Working-directory convention

All scripts read their inputs with **relative paths** (e.g. `read_csv("derf.csv")`,
`stan(file = "Bayesian_contamination_model.stan")`, `readRDS(...)`). Before running a
script, set the R working directory to the folder that contains the files it needs.
The exact file names below include spaces, so keep the quotes when calling `source()`.

## Usage

### 1. Main analysis

```r
setwd("code/data_and_results")
source("Generating results for figures.R")   # fits the model and writes model_fit.rds
```

### 2. Figures

The figure scripts require `derf.csv` and the fitted model object `model_fit.rds` to be
present in the working directory. Copy `derf.csv` and `model_fit.rds` from
`code/data_and_results/` into `code/figures/` (or adjust `setwd()` / the paths), then:

```r
setwd("code/figures")
source("Figure 1.R")   # repeat for "Figure 2.R", ..., "Figure 5.R"
```

Notes on the figure inputs:

- `Figure 2.R` and `Figure 3.R` load the fitted object `model_fit.rds` produced by the
  main analysis script.
- `Figure 1.R` and `Figure 5.R` load saved posterior-sample objects
  (`standard_unknown_4para_model1_mixture_proportion_more_diffuse.rds` and
  `model_fit_newdata_no_initial_value_250.rds`). These are large MCMC fit objects
  generated during the analysis and are not bundled in this archive; regenerate them by
  running the corresponding model fit and saving the object under the expected file name.

### 3. Simulation studies

Each simulation script runs a **single replicate** for demonstration. To reproduce the
full results reported in the manuscript (500 replicates), change
`number_of_iteration <- 1` to `number_of_iteration <- 500` in each script. The simulation
scripts also fit `Bayesian_contamination_model.stan` and `initial_model.stan`, so keep
copies of those `.stan` files in the working directory.

```r
setwd("code/simulation/simulation_1"); source("simulation_1_data_generation.R")
setwd("code/simulation/simulation_2"); source("simulation_2_data_generation.R")
```

## Data

`data/Serial_Dilution_Assay_Sample_Data.csv` provides a sample dataset. The original study
data require IRB approval; please contact the corresponding author for access.
