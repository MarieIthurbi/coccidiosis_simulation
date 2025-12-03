################################################################################
# Coccidiosis SEIS Epidemiological Model - Main Simulation Function
################################################################################
# This script implements a stochastic individual-based model (IBM) combining
# the Gillespie algorithm for disease transitions with deterministic processes
# for growth and environmental contamination.
#
# Model structure: Susceptible -> Exposed -> Infected -> Susceptible (SEIS)
# - Multiple reinfections possible with acquired immunity
# - Individual heterogeneity in susceptibility, shedding, immunity, and growth
# - Environmental transmission through oocyst contamination
# - Growth impacts with tolerance and compensatory growth mechanisms
#
# Author: [Marie Ithurbide _ INRAE]
# Date: [03/15/2025]
# Version: 1.0
################################################################################

################################################################################
# Main Simulation Function: cocci_simu
################################################################################
# Simulates coccidiosis transmission dynamics in a closed poultry population
# using a hybrid stochastic-deterministic framework.
#
# The function requires pre-defined global parameters (loaded via parameter file)
# Returns a list containing multiple output datasets for analysis
################################################################################

cocci_simu <- function() {
  
  # Record simulation timestamp for tracking
  Simulation_date <- Sys.time() |> as.character()

  ## Initialize Data Storage Structures
  ## ---------------------------------------------------------------------------

  SIRdata <- list()      # Disease status counts over time
  SIRdata_env <- list()  # Environmental contamination over time
  PERFs <- list()        # Individual performance metrics (weight, status)
  nReps <- reps          # Number of simulation replicates

  ## Main Simulation Loop
  ## ---------------------------------------------------------------------------
  # Each iteration represents one independent epidemic realization

  for (rep in 1:nReps) {

    # Population initialization
    N <- Ntot              # Total number of animals
    Area <- Area           # Floor area (m²) for density calculations
    pop <- data.table(ID = 1:N)  # Initialize population with unique IDs
    T_list <- c(0)         # Time tracking list

    ## Assign Individual Phenotypes
    ## -------------------------------------------------------------------------
    # Each animal receives random phenotypic values from specified distributions
    # Values are on log-scale (will be exponentiated in calculations)

    pop[, `:=`(
      g = rnorm(.N, mu_g, Sig_g),        # Log-susceptibility
      f = rnorm(.N, mu_f, Sig_f),        # Log-excretion potential
      imm = rnorm(.N, mu_imm, Sig_imm),  # Log-recoverability
      K = rnorm(.N, mu_K, Sig_K),        # Growth potential (actual scale)
      t = rnorm(.N, mu_t, Sig_t),        # Log-tolerance
      r = rnorm(.N, mu_r, Sig_r)         # Log-resilience
    )]

    ## Initialize Body Weights
    ## -------------------------------------------------------------------------
    # Actual weight considers age at exposure and individual growth potential
    # Unperturbed weight tracks what weight would be without infection

    pop[, weight := Birth_weight + Age_exp * K / 80]
    pop[, unperturbed_weight := Birth_weight + Age_exp * K / 80]

    ## Initialize Disease Status and Tracking Variables
    ## -------------------------------------------------------------------------

    # Define health status levels: Susceptible, Exposed, Infected, Dead
    SEI_lvls <- c("S", "E", "I", "D")
    pop[, status := factor("S", SEI_lvls)]  # All start susceptible

    # Initialize infection history tracking variables
    # ninf: number of infections experienced
    # Texp#: time of exposure for infection #
    # Tinf#: time infection # became active (shedding starts)
    # Trec#: time of recovery from infection #
    # Tdeath: time of death (if applicable)
    pop[, `:=`(
      ninf = 0L,
      Texp1 = NA_real_, Texp2 = NA_real_, Texp3 = NA_real_, Texp4 = NA_real_,
      Tinf1 = NA_real_, Tinf2 = NA_real_, Tinf3 = NA_real_, Tinf4 = NA_real_,
      Trec1 = NA_real_, Trec2 = NA_real_, Trec3 = NA_real_, Trec4 = NA_real_,
      Tdeath = NA_real_
    )]

    ## Seed Initial Infection
    ## -------------------------------------------------------------------------
    # Start epidemic with one exposed animal (index case)

    pop[1, `:=`(  
      status = "E",
      Tinf1 = 0,     # Will become infectious at time = incubation period
      ninf = 1L
    )]
    
    # Time step index for data recording
    ts_idx <- 1L

    ## Initialize Population Status Counters
    ## -------------------------------------------------------------------------
    time <- 0
    S <- N - 1L    # Susceptible: all except index case
    E <- 1L        # Exposed: the index case
    I <- 0L        # Infected: none yet (index case in incubation)
    D <- 0L        # Dead: none yet

    # Record initial status counts
    SIRts <- list()
    SIRts[[ts_idx]] <- data.table(time, S, E, I, D)

    ## Initialize Environmental Contamination
    ## -------------------------------------------------------------------------
    nInfectFeces <- 0L  # No oocysts in environment at start

    SIRts_env <- list()
    SIRts_env[[ts_idx]] <- data.table(time, nInfectFeces)


    ############################################################################
    # GILLESPIE ALGORITHM - Stochastic Event-Driven Simulation
    ############################################################################
    # The Gillespie algorithm simulates discrete events (disease transitions)
    # at random times. Between events, deterministic processes (growth,
    # environmental decay) occur continuously.
    #
    # Algorithm steps:
    # 1. Calculate event rates for all individuals
    # 2. Determine time to next event (exponential distribution)
    # 3. Execute event OR deterministic update (whichever comes first)
    # 4. Update system state and repeat
    ############################################################################

    time <- 0
    while (time < Tmax - Age_exp) {
      ts_idx <- ts_idx + 1L

      ## Calculate Individual Event Rates
      ## -----------------------------------------------------------------------
      # Each animal's transition rate depends on their status and phenotype

      # SUSCEPTIBLE → EXPOSED (infection)
      # Rate increases with: environmental contamination, individual susceptibility
      # Rate decreases with: previous infections (immunity c1), area (dilution)
      pop[status == "S", eRates := TRbeta1 * exp(g) * nInfectFeces / ((c1^ninf) * Area)]

      # EXPOSED → INFECTED (end of incubation)
      # Fixed rate (exponentially distributed incubation period)
      pop[status == "E", eRates := RRalpha1]

      # INFECTED → SUSCEPTIBLE (recovery)
      # Rate increases with: recoverability, number of previous infections
      pop[status == "I", eRates := RRalpha3 * exp(imm) * c2^(ninf - 1)]

      # DEAD animals have no transitions
      pop[status == "D", eRates := 0]      
      
      ## Determine Time Step and Event Type
      ## -----------------------------------------------------------------------
      # Compete stochastic events against deterministic time step
      # The process with shortest time occurs first

      eRates_stoch <- sum(pop$eRates)  # Total rate of all stochastic events

      # Calculate time to next event
      dts <- c(
        stoch = if (eRates_stoch > 0) rexp(1L, eRates_stoch) else 1e10,  # Exponential
        det   = dT  # Fixed deterministic time step
      )

      dt <- min(dts, na.rm = TRUE)        # Shortest time wins
      which_dt <- names(which.min(dts))   # Which process occurs?
      time <- time + dt                   # Advance time
      T_list <- c(T_list, time)           # Record time point

      ## Execute Stochastic Event (if selected)
      ## -----------------------------------------------------------------------
      if (which_dt == "stoch") {

        # Select which individual experiences the event
        # Probability proportional to individual event rate
        IDnextT <- pop[, sample(.N, 1, prob = eRates)]
        status <- pop$status[IDnextT]

        # Execute the appropriate transition based on current status
        # ----------------------------------------------------------

        if (status == "S") {
          # INFECTION EVENT: S → E
          set(pop, IDnextT, "status", "E")
          S <- S - 1L
          E <- E + 1L

          # Track infection number and exposure time
          pop[ID == IDnextT, ninf := ninf + 1L]
          pop[ID == IDnextT & ninf == 1, Texp1 := time]
          pop[ID == IDnextT & ninf == 2, Texp2 := time]
          pop[ID == IDnextT & ninf == 3, Texp3 := time]
          pop[ID == IDnextT & ninf == 4, Texp4 := time]
        }

        else if (status == "E") {
          # INCUBATION COMPLETE: E → I (starts shedding oocysts)
          set(pop, IDnextT, "status", "I")
          I <- I + 1L
          E <- E - 1L

          # Record time infection became active
          pop[ID == IDnextT & ninf == 1, Tinf1 := time]
          pop[ID == IDnextT & ninf == 2, Tinf2 := time]
          pop[ID == IDnextT & ninf == 3, Tinf3 := time]
          pop[ID == IDnextT & ninf == 4, Tinf4 := time]
        }

        else if (status == "I") {
          # RECOVERY EVENT: I → S (susceptible to reinfection)
          set(pop, IDnextT, "status", "S")
          I <- I - 1L
          S <- S + 1L

          # Record recovery time
          pop[ID == IDnextT & ninf == 1, Trec1 := time]
          pop[ID == IDnextT & ninf == 2, Trec2 := time]
          pop[ID == IDnextT & ninf == 3, Trec3 := time]
          pop[ID == IDnextT & ninf == 4, Trec4 := time]
        }

      } 

      # Record status counts
      SIRts[[ts_idx]] <- data.table(time, S, E, I, D)

      ############################################################################
      # DETERMINISTIC PROCESSES (occur every time step)
      ############################################################################
      # These continuous processes occur regardless of stochastic events
      ############################################################################

      ## Environmental Oocyst Dynamics
      ## -----------------------------------------------------------------------
      # Change in oocyst load = shedding - decay
      # Shedding: sum of all infected individuals' contributions (modified by f)
      # Decay: exponential degradation at rate RRgamma4

      nInfectFeces <- nInfectFeces - nInfectFeces * RRgamma4 * dt +
        sum(pop[status %in% c("I"), exp(f)]) * Rshed * dt

      ## Body Weight Dynamics
      ## -----------------------------------------------------------------------
      # Weight change depends on health status:
      # - Dead: weight = 0
      # - Infected: reduced growth (tolerance t reduces impact)
      # - Susceptible/Exposed: normal growth + compensatory catch-up

      pop[, weight := fifelse(
        status == "D",
        0,  # Dead animals have zero weight

        fifelse(
          status == "I",
          # Infected animals: reduced growth during infection
          # Growth = (K/80) * [1 - (RRgamma5 / exp(t))] * dt
          # Higher tolerance (t) maintains more growth during infection
          weight + (K / 80) * (1 - (RRgamma5 / exp(t))) * dt,

          # Non-infected animals: normal growth + compensatory growth
          # Compensatory term: exp(r) * c4 * [1 - (actual/expected)]
          # Accelerates growth when weight deficit exists
          weight + (K / 80) * dt * (1 + exp(r) * c4 * (1 - (weight / unperturbed_weight)))
        )
      )]

      # Update expected weight (no infection scenario)
      pop[, unperturbed_weight := unperturbed_weight + (K / 80) * dt]

      ## Mortality Check
      ## -----------------------------------------------------------------------
      # Animals die if weight falls below threshold proportion of expected weight

      IDdead <- pop$ID[(pop$weight / pop$unperturbed_weight) < lim_death & pop$status != "D"]

      if (length(IDdead) > 0) {
        set(pop, IDdead, "status", "D")
        I <- I - 1L  # Assumes dead animals were infected (could be more general)
        D <- D + 1L
        pop[ID %in% IDdead, Tdeath := time]
        pop[ID %in% IDdead, weight := 0]
      }

      ## Record Performance Data
      ## -----------------------------------------------------------------------
      this_time <- time
      pop[, time := this_time]
      PERFs <- rbind(PERFs, pop[, .(rep, ID, time, status, weight, unperturbed_weight)])
      SIRts_env[[ts_idx]] <- data.table(time, nInfectFeces)

    } # End of while loop (simulation complete for this replicate)
    
    ############################################################################
    # Finalize Replicate Data
    ############################################################################

    ts_idx <- ts_idx + 1

    # Record final time point
    SIRts[[ts_idx]] <- data.table(time = Tmax - Age_exp, S, E, I, D)
    SIRts <- rbindlist(SIRts)  # Combine all time steps

    SIRts_env[[ts_idx]] <- data.table(time = Tmax - Age_exp, nInfectFeces)
    SIRts_env <- rbindlist(SIRts_env)

    # Store replicate results
    SIRdata[[rep]] <- list(pop = pop, SIRts = SIRts)
    SIRdata_env[[rep]] <- list(feces = nInfectFeces, SIRts_env = SIRts_env)

  } # End of replicates loop


  ################################################################################
  # POST-PROCESSING: Prepare Output Datasets
  ################################################################################

  ## Format Status Data for Plotting
  ## ---------------------------------------------------------------------------

  # Combine all replicates
  dfplot <- SIRdata |> lapply(\(x) x$SIRts) |> rbindlist(idcol = "rep")

  # Reshape to long format for ggplot
  dfplot2 <- data.table(
    time = rep(dfplot$time, 4),
    rep = rep(dfplot$rep, 4),
    Status = rep(SEI_lvls, each = nrow(dfplot)),
    var = unlist(dfplot[, ..SEI_lvls])
  )

  ## Format Environmental Data
  ## ---------------------------------------------------------------------------

  dfplot2_env <- SIRdata_env |> lapply(\(x) x$SIRts_env) |> rbindlist(idcol = "rep")
  dfplot2_env$Status <- "Feces"

  ## Calculate Weight Deviation Metrics
  ## ---------------------------------------------------------------------------
  # Weight deviation = percentage below expected weight

  PERFs[, weight_deviation := 100 * (1 - weight / unperturbed_weight)]

  # Set dead animals' metrics to NA (not zero) for proper plotting
  PERFs$weight_deviation[PERFs$weight == 0] <- NA
  PERFs$weight[PERFs$weight == 0] <- NA


  ################################################################################
  # INDIVIDUAL-LEVEL SUMMARY DATASET (ANI_SUMMARY)
  ################################################################################
  # Compile detailed infection history and outcomes for each animal
  ################################################################################

  # Extract final population state from all replicates
  ANI_SUMMARY <- SIRdata |> lapply(\(x) x$pop) |> rbindlist(idcol = "rep")

  ## Calculate Duration Metrics
  ## ---------------------------------------------------------------------------

  ANI_SUMMARY[, `:=`(
    # Incubation periods for each infection episode
    exposed_duration_1 = Tinf1 - Texp1,
    exposed_duration_2 = Tinf2 - Texp2,
    exposed_duration_3 = Tinf3 - Texp3,
    exposed_duration_4 = Tinf4 - Texp4,

    # Duration of active infection (shedding period) for each episode
    infection_duration_1 = Trec1 - Tinf1,
    infection_duration_2 = Trec2 - Tinf2,
    infection_duration_3 = Trec3 - Tinf3,
    infection_duration_4 = Trec4 - Tinf4,

    # Time between recovery and next infection (susceptible period)
    recovery_duration_1 = Tinf2 - Trec1,
    recovery_duration_2 = Tinf3 - Trec2,
    recovery_duration_3 = Tinf4 - Trec3,

    # Record simulation timestamp
    Time_Simu = Simulation_date
  )]

  ## Adjust Infection Duration for Animals That Died
  ## ---------------------------------------------------------------------------
  # If animal died during infection, duration is until death (not recovery)

  ANI_SUMMARY[!is.na(Tdeath) & ninf == 1, `:=`(
    infection_duration_1 = Tdeath - Tinf1
  )]

  ANI_SUMMARY[!is.na(Tdeath) & ninf == 2, `:=`(
    infection_duration_2 = Tdeath - Tinf2
  )]

  ANI_SUMMARY[!is.na(Tdeath) & ninf == 3, `:=`(
    infection_duration_3 = Tdeath - Tinf3
  )]

  ANI_SUMMARY[!is.na(Tdeath) & ninf == 4, `:=`(
    infection_duration_4 = Tdeath - Tinf4
  )]
  
  ## Calculate Daily Average Performance Metrics
  ## ---------------------------------------------------------------------------

  PERFs[, Day := round(time, 0)]
  PERFs2 <- PERFs[, lapply(.SD, mean, na.rm = TRUE),
                  .SDcols = c("weight", "unperturbed_weight", "weight_deviation"),
                  by = .(ID, rep, Day, status)]

  ## Add Weight Deviation Summary to Individual Data
  ## ---------------------------------------------------------------------------

  ANI_SUMMARY[, `:=`(
    max_weight_deviation = PERFs2[, max(weight_deviation, na.rm = TRUE), .(ID, rep)][, V1],
    sum_weight_deviation = PERFs2[, sum(weight_deviation, na.rm = TRUE), .(ID, rep)][, V1]
  )]

  ## Calculate Total Time Infected
  ## ---------------------------------------------------------------------------

  ANI_SUMMARY$d_inf <- rowSums(ANI_SUMMARY[, c("infection_duration_1",
                                               "infection_duration_2",
                                               "infection_duration_3",
                                               "infection_duration_4")],
                               na.rm = TRUE)

  ## Mortality Indicator
  ## ---------------------------------------------------------------------------

  ANI_SUMMARY$death <- 0
  ANI_SUMMARY$death[ANI_SUMMARY$weight == 0] <- 1


  ################################################################################
  # SIMULATION PARAMETER RECORD (SIMULATION_param)
  ################################################################################
  # Store all simulation settings for reproducibility and metadata
  ################################################################################

  SIMULATION_param <- data.table(
    Time_Simu = Simulation_date,
    Ntot = Ntot,
    reps = reps,
    Tmax = Tmax,
    Age_exp = Age_exp,
    mu_g = mu_g,
    Sig_g = Sig_g,
    mu_f = mu_f,
    Sig_f = Sig_f,
    mu_imm = mu_imm,
    Sig_imm = Sig_imm,
    mu_K = mu_K,
    Sig_K = Sig_K,
    mu_t = mu_t,
    Sig_t = Sig_t,
    mu_r = mu_r,
    Sig_r = Sig_r
  )
  
  ################################################################################
  # TIME SERIES SUMMARY DATA (dfplot2_env with prevalence and weight metrics)
  ################################################################################

  dfplot2_env[, Day := round(time, 0)]

  ## Calculate Prevalence (proportion of living animals infected)
  ## ---------------------------------------------------------------------------
  dfplot2_env$prev <- dfplot2$var[dfplot2$Status == "I"] /
                      (Ntot - dfplot2$var[dfplot2$Status == "D"])

  ## Add Weight Deviation Time Series
  ## ---------------------------------------------------------------------------
  # Mean and SD of weight deviation at each time point
  AA <- cbind(
    aggregate(weight_deviation ~ time + rep, data = PERFs, mean),
    aggregate(weight_deviation ~ time + rep, data = PERFs, sd)[, 3]
  )
  colnames(AA)[3:4] <- c("mean_wd", "sd_wd")
  dfplot2_env <- left_join(dfplot2_env, AA)


  ################################################################################
  # REPLICATE-LEVEL SUMMARY DATASET (SIMULATION)
  ################################################################################
  # Aggregate statistics for each simulation replicate
  ################################################################################

  ## Define Variables to Average Across Individuals
  ## ---------------------------------------------------------------------------
  cols <- c("ninf", "unperturbed_weight", "max_weight_deviation",
            "sum_weight_deviation", "d_inf", "Texp1", "exposed_duration_1",
            "exposed_duration_2", "exposed_duration_3", "exposed_duration_4",
            "infection_duration_1", "infection_duration_2", "infection_duration_3",
            "infection_duration_4", "recovery_duration_1", "recovery_duration_2",
            "recovery_duration_3")

  SIMULATION <- ANI_SUMMARY[, lapply(.SD, mean, na.rm = TRUE),
                            .SDcols = cols,
                            by = .(Time_Simu, rep)]

  ## Add Final Weight Statistics
  ## ---------------------------------------------------------------------------
  # Mean weight at slaughter (excluding dead animals)
  weight_means <- ANI_SUMMARY[weight > 0, .(weight = mean(weight, na.rm = TRUE)),
                              by = .(Time_Simu, rep)]

  SIMULATION <- merge(SIMULATION, weight_means, by = c("Time_Simu", "rep"))

  # Standard deviation of final weights
  SIMULATION$final_weight_sd <- aggregate(weight ~ Time_Simu + rep,
                                          data = subset(ANI_SUMMARY, weight > 0),
                                          sd)[, 3]

  ## Add Mortality Statistics
  ## ---------------------------------------------------------------------------
  SIMULATION$n_death <- aggregate(death ~ Time_Simu + rep,
                                  data = ANI_SUMMARY, sum)[, 3]
  SIMULATION$prop_death <- SIMULATION$n_death /
                           rep(SIMULATION_param$Ntot, each = SIMULATION_param$reps)

  ## Add Environmental and Economic Metrics
  ## ---------------------------------------------------------------------------
  SIMULATION[, `:=`(
    max_feces_env = dfplot2_env[, max(nInfectFeces, na.rm = TRUE), by = rep][, V1],
    mean_feces_env = dfplot2_env[, mean(nInfectFeces, na.rm = TRUE), by = rep][, V1],
    mean_prev = dfplot2_env[, mean(prev, na.rm = TRUE), by = rep][, V1],
    # Meat production loss (proportion of potential production lost)
    prop_final_meet_loss = 1 - (ANI_SUMMARY[, sum(weight, na.rm = TRUE), by = rep][, V1] /
                                ANI_SUMMARY[, sum(unperturbed_weight, na.rm = TRUE), by = rep][, V1])
  )]


  ################################################################################
  # RETURN OUTPUT LIST
  ################################################################################
  # Package all results for user analysis and visualization
  ################################################################################

  return(list(
    SIMULATION = SIMULATION,                # Replicate-level summaries
    SIMULATION_param = SIMULATION_param,    # Simulation parameters
    PERFs = PERFs,                          # Individual performance time series
    ANI_SUMMARY = ANI_SUMMARY,              # Individual infection histories
    dfplot2 = dfplot2,                      # Status counts time series (long format)
    dfplot2_env = dfplot2_env               # Environmental and prevalence time series
  ))

} # End of cocci_simu function

################################################################################
# END OF FILE
################################################################################
