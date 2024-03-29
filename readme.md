Mortality risk from United States coal electricity generation
================
Lucas Henneman
2023-08-11

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

This repository houses code for the manuscript entitled “Mortality risk
from United States coal electricity generation”.

# Directory contents

-   `code` This directory stores all code used as part of the project.
    -   `00_organize_data.R` Code in this file reads in the outcome and
        exposure data into consistent formats.
    -   `01_coal_pm_mortality_models.R` Code in this file runs the
        stratified poisson models.
    -   `02_CI_bootstrap.R` Code in this file calculates the
        bootstrapped confidence intervals on the estimated hazard ratios
    -   `03_deaths_contributed.R` Code in this file combines coal power
        plant information with poisson model output to calculate deaths
        associated with each coal power plant’s SO<sub>2</sub>
        emissions.
    -   `04_plots.R` Code in this file creates many of the charts in
        figures in the manuscript and supplementary text.
    -   `05_compare_geos_chem.R` Code in this file compares the power
        plant deaths calculated in this study with comparable deaths
        calculated using GEOS-Chem adjoint sensitivities
    -   `06_change_over_change.R` Code in this file conducts the
        first-differences analysis described in the supplementary text
    -   `07_plot_hyads_biggest_impactors.R` Code in this file creates
        figures presented in the supplementary text showing facilities
        that contribute to coal PM<sub>2.5</sub> in metropolitan
        statistical areas
-   `data` Sensitive health data is not included in this repository, but
    all other datasets are saved here.
    -   `cache_data` This sub-directory houses various ZIP code level
        data for input into the models. Sensitive data and large files
        have not been pushed to Github.
    -   `adjoint_results` This sub-directory houses deaths (total
        population) attributable to each coal unit from the GEOS-Chem
        adjoint sensitivities modeling.
    -   `coalpm25_to_obs_adjustment.csv` This file contains regional and
        year-specific adjustment factors for HyADS-derived coal
        PM<sub>2.5</sub> to regional observed sulfate concentrations
-   `figures` This directory contains figures from the manuscript
-   `results` This directory contains two files related to the Risk
    Ratio estimations in the manuscript
    -   `poisson_model_coefs.csv` Contains the central estimates of the
        relevant coefficients from the Poisson models used to estimate
        the RR
    -   `loglinear_coefs_boots.RData` Contains the bootstrapped
        estimates of the relevant coefficients from the Poisson models
        used to estimate the RR confidence intervals
