################################################################################
# Model Parameters for Coccidiosis SEIS Model
################################################################################
# This file defines all epidemiological and biological parameters for the
# coccidiosis transmission model in poultry. The model simulates infection
# dynamics, oocyst shedding, immune response, and growth impacts.
################################################################################

## Epidemiological Parameters
## ---------------------------------------------------------------------------

# Transmission rate parameter (r1)
# Controls the rate of transition from Susceptible (S) to Exposed (E)
# Higher values indicate increased susceptibility to infection
TRbeta1 <- 1.61

# Oocyst shedding rate (r4)
# Number of oocysts shed per day per infected animal
# This determines the environmental contamination level
Rshed <- 2.84

# Environmental decay rate (r3)
# Rate at which oocysts degrade in the environment
# Higher values mean faster decay of infectious oocysts
RRgamma4 <- 1.15

# Incubation period parameter (d)
# Inverse of the duration from Exposed (E) to Infected (I)
# RRalpha1 = 1/4 means an average 4-day incubation period
RRalpha1 <- 1/4

# Immunity effect on susceptibility (c1)
# Multiplicative factor reducing reinfection risk
# Higher values indicate stronger protection from previous infections
c1 <- 5.65

## Growth and Physiological Parameters
## ---------------------------------------------------------------------------

# Initial body weight at start of study (grams)
Birth_weight <- 140

# Time step for deterministic processes (days)
# Smaller values increase precision but computational cost
dT <- 0.1

# Growth reduction coefficient (c3)
# Effect of infection on daily weight gain, strain-specific
# Lower values indicate greater growth impact
RRgamma5 <- 0.4

# Recovery rate parameter (r2)
# Inverse of mean infection duration for first infection
# RRalpha3 = 1/12.1 means average 12.1 days of first infection
RRalpha3 <- 1/12.1

# Immunity effect on infection duration (c2)
# Multiplicative factor reducing duration of subsequent infections
# Higher values indicate faster recovery with acquired immunity
c2 <- 2.3

# Compensatory growth rate (c4)
# Speed at which animals recover lost weight after infection
# Higher values indicate faster catch-up growth
c4 <- 0.1

# Mortality threshold
# Proportion of expected weight below which death occurs
# lim_death = 0.8 means death when actual weight < 80% of expected
lim_death <- 0.8

## Individual Heterogeneity Parameters
## ---------------------------------------------------------------------------
# These parameters define the distribution of phenotypic variation among
# individuals. Each trait follows a normal distribution with mean (mu) and
# standard deviation (Sig). Setting Sig = 0 creates a homogeneous population.

# Susceptibility heterogeneity (g)
# Individual variation in infection risk
# Higher g increases individual susceptibility
mu_g <- 0      # Mean log-susceptibility
Sig_g <- 0     # SD of log-susceptibility

# Excretion potential heterogeneity (f)
# Individual variation in oocyst shedding capacity
# Higher f increases individual shedding rate
mu_f <- 0      # Mean log-excretion
Sig_f <- 0     # SD of log-excretion

# Recoverability heterogeneity (imm)
# Individual variation in infection duration
# Higher imm accelerates recovery 
mu_imm <- 0    # Mean log-recoverability
Sig_imm <- 0   # SD of log-recoverability

# Growth potential heterogeneity (K)
# Individual variation in target body weight at 80 days (grams)
# Determines maximum achievable body weight
mu_K <- 1800   # Mean final weight potential (g)
Sig_K <- 0     # SD of final weight potential

# Disease tolerance heterogeneity (t)
# Individual variation in ability to maintain growth during infection
# Higher t reduces the negative impact of infection on weight gain
mu_t <- 0      # Mean log-tolerance
Sig_t <- 0     # SD of log-tolerance

# Disease resilience heterogeneity (r)
# Individual variation in compensatory growth capacity
# Higher r increases catch-up growth rate after infection
mu_r <- 0      # Mean log-resilience
Sig_r <- 0     # SD of log-resilience













