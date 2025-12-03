# Coccidiosis SEIS Model

**Author:** Marie Ithurbide (INRAE)
**Version:** 1.0
**Date:** March 2025

## Overview

Stochastic individual-based model (IBM) for coccidiosis transmission in poultry. Combines Gillespie algorithm for disease transitions with deterministic processes for growth and environmental contamination.

**Model:** Susceptible → Exposed → Infected → Susceptible (SEIS)
- Multiple reinfections with acquired immunity
- Individual heterogeneity (susceptibility, shedding, recoverability, tolerance, resilience)
- Environmental transmission via oocysts
- Growth impacts with compensatory mechanisms

---

## Files

### 1. `Model.parameters - simple-V1.R`
Defines all model parameters:
- Transmission rates and shedding dynamics
- Incubation and recovery periods
- Immunity effects
- Growth parameters and mortality threshold
- Individual heterogeneity distributions

### 2. `Function Modele coccidiose_simple V1.R`
Main simulation engine with `cocci_simu()` function.

**Process:**
- Initializes population with individual phenotypes
- Seeds initial infection
- Runs Gillespie algorithm for stochastic transitions (S→E→I→S)
- Updates deterministic processes: environmental dynamics, growth, mortality
- Tracks infection history (up to 4 infections per animal)

**Returns 6 datasets:**
- `SIMULATION`: Replicate summaries (weights, mortality, prevalence)
- `SIMULATION_param`: Parameter record
- `PERFs`: Individual time series
- `ANI_SUMMARY`: Individual infection histories
- `dfplot2`: Disease status counts
- `dfplot2_env`: Environmental and prevalence time series

### 3. `Coccidiose - Run the model - simple.Rmd`
R Markdown demonstration file with:
- Package loading and parameter configuration
- Simulation execution
- Visualizations (epidemic dynamics, growth curves, relationships)
- Summary statistics

**Usage:** Open in RStudio, modify parameters, knit to HTML.

---

## Quick Start

```r
# 1. Load packages
library(data.table)
library(ggplot2)
library(dplyr)

# 2. Source files
source('Model.parameters - simple-V1.R')
source('Function Modele coccidiose_simple V1.R')

# 3. Configure simulation
reps <- 10
Ntot <- 20
Tmax <- 84
Age_exp <- 10
Area <- 20

# 4. Set heterogeneity (0 = homogeneous)
mu_g <- 0; Sig_g <- 0       # Susceptibility
mu_f <- 0; Sig_f <- 0       # Excretion
mu_imm <- 0; Sig_imm <- 0   # Recoverability
mu_K <- 1800; Sig_K <- 0    # Growth potential
mu_t <- 0; Sig_t <- 0       # Tolerance
mu_r <- 0; Sig_r <- 0       # Resilience

# 5. Run
results <- cocci_simu()

# 6. View results
summary(results$SIMULATION$weight)
summary(results$SIMULATION$prop_death)
```

---

## Requirements

**R:** ≥ 4.0.0

**Packages:**
- Required: `data.table`, `ggplot2`, `dplyr`
- Optional (for Rmd): `ggpubr`, `viridis`, `rockchalk`, `gt`

```r
install.packages(c("data.table", "ggplot2", "dplyr"))
```

---

## Contact

Marie Ithurbide - INRAE
