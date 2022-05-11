Coal transition and health
================
Lucas Henneman
05/10/2022

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

This repository houses code for the manuscript entitled “Coal’s Toll: 22
years of Medicare deaths attributable to electricity generation”.

# Directories

-   `code` This directory stores all code used as part of the project.
    -   `00_organize_data.R` Code in this file reads in the outcome and
        exposure data into consistent formats.
    -   `01_coal_pm_mortality_models.R` Code in this file runs the
        stratefied poisson models.
    -   `02_CI_bootstrap.R` Code in this file calculates the
        bootstrapped confidence intervals on the estimated hazard ratios
    -   `03_deaths_contributed.R` Code in this file combines coal power
        plant information with poisson model output to calculate deaths
        associated with each coal power plant’s SO<sub>2\<\> emissions.
    -   `04_plots.R` Code in this file creates many of the charts in
        figures in the manuscript.
    -   `05_compare_geos_chem.R`Code in this file compares the power
        plant deaths calculated in this study with comperable deaths
        calculated using GEOS-Chem adjoint sensitivities
-   `data` Sensitive health data is not included in this repository, but
    all other datasets are saved here.
    -   `admissions` This subdirectory houses raw admissions data
    -   `cache_data` This subdirectory houses monthly and annual ZIP
        code rates of multiple health outcomes
    -   `denom` This subdirectory houses monthly denominator files for
        all Medicare beneficiaries
    -   `denom_annual` This subdirectory houses annual denominator files
        for all Medicare beneficiaries
    -   `denom_ffs` This subdirectory houses monthly denominator files
        for all fee for service beneficiaries
-   `figures` This directory contains figures from the manuscript
-   `results` This directory contains two files
