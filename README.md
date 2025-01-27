# **Modelling the Spatial Dependence of Multi-species Point Patterns**

This repository contains the implementation and analysis associated with the study comparing two multivariate point process models: 
- **Multivariate Log Gaussian Cox Process (MLGCP)**
-  **Saturated Pairwise Interaction Gibbs Point Process (SPIGPP)**

The study evaluates their efficacy in modeling spatial point patterns in ecology, focusing on their ability to handle clustering and regulation under varying conditions.

---

## **Overview**

The spatial analysis of ecological data, such as the observed locations of trees, nests, or animal sightings, is crucial in understanding environmental patterns. This study compares MLGCP and SPIGPP, highlighting their strengths and weaknesses when prior knowledge of driving mechanisms is limited. Using synthetic and real datasets, we examine their predictive accuracy for the empirical K function and evaluate their performance under different interaction types and intensities.

- **MLGCP**: Effective for accounting for complex, unobserved spatial heterogeneities.
- **SPIGPP**: Better at estimating interactions (even in misspecified cases).

For detailed results, please refer to the [pre-print](https://doi.org/10.22541/au.172623817.72677555/v1) 📄.

## **Repository Contents**
- `code/`: Updated code for simulating from a MLGCP.
- `src/`: Code to fit MLGCP and SPIGPP models and evaluate their performance.
- `README.md`: Overview and instructions for using the repository.


## **Installation**

Ensure you have R (version 4.3.1 or above) installed on your system. Required R packages include those available in CRAN, along with the Multilogreg package (see below for instructions).

1. Install Required R Packages:
To install the required packages, run the provided script: [Install required packages.R](Install%20required%20packages.R)

2. The original [Multilogreg package](https://github.com/kristianhessellund/Multilogreg.git) has compatibility and dependency issues. Use the fixed version in [IbTJensen](https://github.com/IbTJensen/Multilogreg.git).

## **Usage**

To replicate the analysis:
1. run the script [sim_MLGCP_adjusted.R](Sim_MLGCP_adjusteed.R). This will enable the use of any given window and multiple covariates when simulating from MLGCP.
2. Try the [reproducible_MLGCP_Scenarios.R](reproducible_MLGCP_Scenarios.R), then,
3. [reproducible_PPJSDM scenarios.R](reproducible_PPJSDM%20scenarios.R)
4. Next try the [reproducible_5_species_simulation.R](reproducible_5_species_simulation.R)
5. Finally try the real data example on Swamp data Analysis [Reproducible Swamp Data Analysis.R](Reproducible%20Swamp%20Data%20Analysis.R)

Results include: Predicted K functions for synthetic and real datasets and model performance metrics for clustering and regulation.

## Acknowledgments

Original Multilogreg Package: Kristian Hessellund's [Multilogreg package](https://github.com/kristianhessellund/Multilogreg.git) 

Fixed Multilogreg Package: [IbT Jensen's version](https://github.com/IbTJensen/Multilogreg.git)

## Pre-print
Chathuri Samarasekara, Yan Wang, Ian Flint. A Comparison of Multivariate Log Gaussian Cox Process and Saturated Pairwise Interaction Gibbs Point Process. Authorea. September 13, 2024. DOI: 10.22541/au.172623817.72677555/v1

For any questions or feedback, feel free to reach out:

Email: [chathuri.samarasekara@rmit.edu.au](mailto:chathuri.samarasekara@rmit.edu.au)

GitHub: [Chathuri Samarasekara](https://github.com/chathuri-sam)



