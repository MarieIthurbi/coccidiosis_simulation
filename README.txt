================================================================================
COCCIDIOSIS SEIS EPIDEMIOLOGICAL MODEL - R CODE
================================================================================

Author: Marie Ithurbide (INRAE)
Date: March 2025
Version: 1.0

================================================================================
OVERVIEW
================================================================================

Stochastic individual-based model (IBM) for simulating coccidiosis transmission
in poultry. Combines Gillespie algorithm for disease transitions with
deterministic processes for growth and environmental contamination.

Model: Susceptible → Exposed → Infected → Susceptible (SEIS)
- Multiple reinfections with acquired immunity
- Individual heterogeneity (susceptibility, shedding, recoverability,
  tolerance, resilience)
- Environmental transmission via oocysts
- Growth impacts with compensatory mechanisms

================================================================================
FILES DESCRIPTION
================================================================================

1. Model.parameters - simple-V1.R
-----------------------------------

Defines all model parameters:
- Transmission rates and shedding dynamics
- Incubation and recovery periods
- Immunity effects
- Growth parameters and mortality threshold
- Individual heterogeneity distributions

Parameters are loaded into global environment for use by the simulation.


2. Function Modele coccidiose_simple V1.R
------------------------------------------

Main simulation engine containing the cocci_simu() function.

Process:
- Initializes population with individual phenotypes
- Seeds initial infection (one exposed animal)
- Runs Gillespie algorithm for stochastic disease transitions (S→E→I→S)
- Updates deterministic processes at each time step:
  * Environmental oocyst dynamics (shedding and decay)
  * Body weight changes (reduced growth during infection, compensatory
    growth after)
  * Mortality based on weight loss
- Tracks infection history for each animal (up to 4 infections)
- Aggregates results across replicates

Requires global variables: reps, Ntot, Tmax, Age_exp, Area, mu_*, Sig_*

Returns 6 datasets:
- SIMULATION: Replicate-level summaries (weights, mortality, prevalence)
- SIMULATION_param: Parameter record
- PERFs: Individual performance time series
- ANI_SUMMARY: Individual infection histories
- dfplot2: Disease status counts (long format)
- dfplot2_env: Environmental and prevalence time series


3. Coccidiose - Run the model - simple.Rmd
-------------------------------------------

R Markdown demonstration file.

Contents:
- Loads required packages
- Sources parameter and function files
- Configures simulation settings
- Runs simulation via cocci_simu()
- Generates visualizations:
  * Multi-panel plots by replicate
  * Summary plots (epidemic dynamics, growth curves)
  * Relationship plots (prevalence vs contamination)
- Produces summary statistics tables

Usage: Open in RStudio, modify parameters if needed, then knit to HTML.

================================================================================
QUICK START
================================================================================

# 1. Load packages
library(data.table)
library(ggplot2)
library(dplyr)

# 2. Source files
source('Model.parameters - simple-V1.R')
source('Function Modele coccidiose_simple V1.R')

# 3. Configure simulation
reps <- 10              # Number of replicates
Ntot <- 20              # Population size
Tmax <- 84              # Final age (days)
Age_exp <- 10           # Age at infection (days)
Area <- 20              # Floor area (m²)

# 4. Set heterogeneity (0 = homogeneous)
mu_g <- 0; Sig_g <- 0       # Susceptibility
mu_f <- 0; Sig_f <- 0       # Excretion
mu_imm <- 0; Sig_imm <- 0   # Recoverability
mu_K <- 1800; Sig_K <- 0    # Growth potential
mu_t <- 0; Sig_t <- 0       # Tolerance
mu_r <- 0; Sig_r <- 0       # Resilience

# 5. Run simulation
results <- cocci_simu()

# 6. View results
summary(results$SIMULATION$weight)       # Final weights
summary(results$SIMULATION$prop_death)   # Mortality

================================================================================
REQUIREMENTS
================================================================================

R >= 4.0.0

Required: data.table, ggplot2, dplyr
Optional (for Rmd): ggpubr, viridis, rockchalk, gt

Install: install.packages(c("data.table", "ggplot2", "dplyr"))

================================================================================
CONTACT
================================================================================

Marie Ithurbide - INRAE

================================================================================
