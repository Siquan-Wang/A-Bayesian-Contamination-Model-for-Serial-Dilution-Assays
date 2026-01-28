# A Bayesian Contamination Model for Serial Dilution Assays

Code for reproducing the results in the manuscript submitted to *Biometrics*.

## Requirements

- R (version 4.0+)
- RStan
- R packages: tidyverse, ggpubr, psych, loo, shinystan, cowplot, Rmpfr, ggrepel

## Usage

### Main Analysis and Figures

```r
source("code/data_and_results/Generating_results_for_figures.R")
source("code/figures/Figure_1.R")  # Repeat for Figures 2-5
```

### Simulation Studies

The simulation code runs a single replicate for demonstration. To reproduce full results (500 replicates), change `number_of_iteration <- 1` to `number_of_iteration <- 500` in each script.

```r
source("code/simulation/simulation_1/simulation_1_data_generation.R")
source("code/simulation/simulation_2/simulation_2_data_generation.R")
```

## Data

Sample data is provided in `data/`. Original study data requires IRB approval; contact the corresponding author for access.
