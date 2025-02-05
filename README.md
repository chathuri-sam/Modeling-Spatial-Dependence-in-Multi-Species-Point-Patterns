# **Modeling Spatial Dependence in Multi-Species Point Patterns**

This repository contains the implementation and analysis associated with the study comparing two multivariate point process models: 
- **Multivariate Log Gaussian Cox Process (MLGCP)**
-  **Saturated Pairwise Interaction Gibbs Point Process (SPIGPP)**

The study evaluates their efficacy in modeling spatial point patterns in ecology, focusing on how well they capture clustering and interaction effects under different conditions.

---

## **Overview**

The spatial analysis of ecological data, such as the observed locations of trees, nests, or animal sightings, is crucial in understanding environmental patterns. This study compares MLGCP and SPIGPP, highlighting their respective strengths and limitations, particularly in cases where the underlying drivers of spatial patterns are unknown. Using synthetic and real datasets, we examine their predictive accuracy assessed using the empirical K function to quantify spatial dependence and evaluate their performance under varying interaction structures.

Summary:

- **MLGCP**: Effective for accounting for complex, unobserved spatial heterogeneities.
- **SPIGPP**: More effective at detecting interactions, even when the true underlying model is misspecified.

For detailed results, please refer to the [pre-print](https://doi.org/10.22541/au.172623817.72677555/v1).

## **Repository Contents**
- `src/`: Code for fitting MLGCP and SPIGPP models in various scenarios, evaluating their performance, and analyzing real data to compare the two models.
- `README.md`: Overview and instructions for using the repository.

---

## **Installation**

1. Ensure that R is installed on your system. The code has been tested on R 4.4.2 and confirmed to work on Linux, Windows, and macOS.
2. Required R packages include those available in CRAN; 
 - *spatstat* - version 3.0-6
 - *ggplot2* - version 3.5.1
 - *ecespa* - version 1.1-17

3. Packages available on Github
  - *RandomFieldsUtils* (fixed) - version 1.2.5
  - *RandomFields* (fixed) - version 3.3.14
  - *PPJSDM* - version 1.0
  - *Multilogreg* (fixed) - version 0.1.0

Note : The original [Multilogreg package](https://github.com/kristianhessellund/Multilogreg.git) does not work since a recent `spatstat` update. Therefore, we use the fixed version in [Multilogreg_updated](https://github.com/chathuri-sam/Multilogreg_updated.git). 
 
To install all required packages, run the provided script: [Install required packages.R](src/Install%20required%20packages.R)

## **Usage**

Follow these steps to reproduce the analysis:
1. Simulate bi-variate MLGCP scenarios and fit them using the SPIGPP models. [reproducible_MLGCP_Scenarios.R](src/reproducible_MLGCP_Scenarios.R), 
2. Next, simulate bi-variate SPIGPP scenarios and fit them using the MLGCP models. [reproducible_PPJSDM scenarios.R](src/reproducible_PPJSDM%20scenarios.R)
   For both of these mis-specified scenarios, we check the model fit by comparing the empirical and fitted K functions.
3. Next, reproduce the five-species simulation study from [Hessellund et al. (2022)](https://doi.org/10.1111/rssc.12530) and fit it with MLGCP and SPIGPP models. Then evaluate the fit using the empirical an fitted K functions [reproducible_5_species_simulation.R](src/reproducible_5_species_simulation.R)
4. Finally, apply both models to the South Carolina Savannah River site dataset from Good & Whipple (1982) using [Reproducible Swamp Data Analysis.R](src/Reproducible%20Swamp%20Data%20Analysis.R).

Detailed results and interpretations are available in the paper.

## Acknowledgments

Original Multilogreg Package: Kristian Hessellund's [Multilogreg package](https://github.com/kristianhessellund/Multilogreg.git) 

## References
1. Good, B. and Whipple, S. (1982). Tree spatial patterns: South Carolina bottomland and swamp forests. Bulletin of the Torrey Botanical Club, 109:529–536.528.
2. Hessellund, K. B., Xu, G., Guan, Y., and Waagepetersen, R. (2022a). Second-Order Semi-Parametric Inference for Multivariate Log Gaussian Cox Processes. Journal of the Royal Statistical Society Series C: Applied Statistics, 71(1):244–268.531. [https://doi.org/10.1111/rssc.12530](https://doi.org/10.1111/rssc.12530)

## Pre-print
Chathuri Samarasekara, Yan Wang, Ian Flint. A Comparison of Multivariate Log Gaussian Cox Process and Saturated Pairwise Interaction Gibbs Point Process. Authorea. September 13, 2024. DOI: 10.22541/au.172623817.72677555/v1

For any questions or feedback, feel free to reach out:

Email: [chathuri.samarasekara@rmit.edu.au](mailto:chathuri.samarasekara@rmit.edu.au)

GitHub: [Chathuri Samarasekara](https://github.com/chathuri-sam)
