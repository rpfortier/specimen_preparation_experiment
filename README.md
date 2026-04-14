# Specimen preservation techniques affect herbarium-derived leaf traits, reflectance spectra, and predictive models 

Contact: Riley P. Fortier (fortier.riley@gmail.com)

This repository contains code and associated data to analyze how drying and alcohol preservation affect leaf functional traits, reflectance spectra, and predictive models built from reflectance spectra. A manuscript is currently in review.

All R scripts run successfully using R software version 4.5.0 on Windows 11 and macOS Tahoe and should work on all standard operating systems. Details of each script including runtimes are included below.

## *Leaf traits and reflectance spectra data*
* specimen_prep_leaf_data.csv
  + Fresh and dry leaf trait data.
* TRY_trait_range.csv
  + Global plant trait data summarized from the TRY plant trait database.
* trees_collected.csv
  + Tree metadata including voucher numbers and species identities.
* dry_hyperspectral_data.csv (download at https://drive.google.com/file/d/1RB4lgGDJq6RVxMK8nTmfETZSyqNUEiSe/view?usp=drive_link)
  + Reflectance values across the entire spectrum (350 nm - 2500 nm) for the abaxial and adaxial side of each leaf.

## *R scripts*
* spectral_functions.R
  + Creates essential functions including %RMSE calculations and derivative calculations
* Spec_prep_trait_analysis.R
  + Cleans data and computes standardized trait values. Analyzes preservation treatment effects on leaf traits and reflectance spectra.
  + This script should take <10 minutes to complete.
* Spec_prep_PLSR.R
  + Code to build partial least-squares regressions for trait prediction.
  + This script takes longer, up to a few hours to complete all models.
* Spec_prep_PLS-DA.R
  + Code to build partial least-squares discriminate analysis for species classification.
  + This script is computationally heavy and requires parallel processing. On Windows 11 with 64 GB of RAM and running 26 cores simultaneously, each PLS-DA model takes several hours to complete. Therefore, multiple days are needed to complete all PLS-DA models.  
