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

For detailed results, please refer to the [pre-print](https://doi.org/10.22541/au.172623817.72677555/v1).

## **Repository Contents**
- `code/`: Updated code for simulating from a MLGCP.
- `src/`: Code to fit MLGCP and SPIGPP models and evaluate their performance.
- `README.md`: Overview and instructions for using the repository.

---

## **Installation**

1. Ensure you have R (version 4.3.1 or above) installed on your system.
2. Required R packages include those available in CRAN; 
 - *spatstat* - version 3.0-6
 - *ggplot2* - version 3.5.1
 - *ecespa* - version 1.1-17
 
3. Packages that are not currently avaiable on CRAN;
  - *RandomFieldsUtils* - version 1.2.5
  - *RandomFields* - version 3.3.14

4. Packages available on Github
  - *PPJSDM* - version 1.0
  - *Multilogreg* - version 0.1.0

Note : The original [Multilogreg package](https://github.com/kristianhessellund/Multilogreg.git) has compatibility and dependency issues. Therefore, use the fixed version in [IbTJensen](https://github.com/IbTJensen/Multilogreg.git).
 
To install all required packages, run the provided script: [Install required packages.R](Install%20required%20packages.R)

## **Usage**

To replicate the analysis:
1. run the script [sim_MLGCP_adjusted.R](Sim_MLGCP_adjusteed.R). This will enable the use of any given window and multiple covariates when simulating from MLGCP (The original version of sim_MLGCP in Multilogreg package uses a square window and a single covariate).
2. Then, simulate bi-variate MLGCP scenarios and fit them using the SPIGPP models. [reproducible_MLGCP_Scenarios.R](reproducible_MLGCP_Scenarios.R), 
3. Next, simulate bi-variate SPIGPP scenarios and fit them using the MLGCP models. [reproducible_PPJSDM scenarios.R](reproducible_PPJSDM%20scenarios.R)
   For both of these mis-specified sceanrios, we check the model fit using the empirical and fitted K functions.
5. Then, replicate the 5 species simulation study in [Hessellund et al. (2022)](https://doi.org/10.1111/rssc.12530) and fit it with SPIGPP models and evaluate the fit using the empirical an fitted K functions [reproducible_5_species_simulation.R](reproducible_5_species_simulation.R)
6. Finally apply MLGCP and SPIGPP models to the South Carolina Savannah river site study (Good & Whipple (1982)) [Reproducible Swamp Data Analysis.R](Reproducible%20Swamp%20Data%20Analysis.R).

The final results and interpretations of them can be found the paper.

## Acknowledgments

Original Multilogreg Package: Kristian Hessellund's [Multilogreg package](https://github.com/kristianhessellund/Multilogreg.git) 

Fixed Multilogreg Package: [IbT Jensen's version](https://github.com/IbTJensen/Multilogreg.git)

## References
1. Good, B. and Whipple, S. (1982). Tree spatial patterns: South Carolina bottomland and swamp forests. Bulletin of the Torrey Botanical Club, 109:529–536.528.
2. Hessellund, K. B., Xu, G., Guan, Y., and Waagepetersen, R. (2022a). Second-Order Semi-Parametric Inference for Multivariate Log Gaussian Cox Processes. Journal of the Royal Statistical Society Series C: Applied Statistics, 71(1):244–268.531. [https://doi.org/10.1111/rssc.12530](https://doi.org/10.1111/rssc.12530)

## Pre-print
Chathuri Samarasekara, Yan Wang, Ian Flint. A Comparison of Multivariate Log Gaussian Cox Process and Saturated Pairwise Interaction Gibbs Point Process. Authorea. September 13, 2024. DOI: 10.22541/au.172623817.72677555/v1

For any questions or feedback, feel free to reach out:

Email: [chathuri.samarasekara@rmit.edu.au](mailto:chathuri.samarasekara@rmit.edu.au)

GitHub: [Chathuri Samarasekara](https://github.com/chathuri-sam)
