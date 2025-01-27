# Modelling the Spatial Dependence of Multi-species Point Patterns

This repository contains the implementation and analysis associated with the study comparing two multivariate point process models: the Multivariate Log Gaussian Cox Process (MLGCP) and the Saturated Pairwise Interaction Gibbs Point Process (SPIGPP). The study evaluates their efficacy in modeling spatial point patterns in ecology, focusing on their ability to handle clustering and regulation under varying conditions.

Overview

The spatial analysis of ecological data, such as the observed locations of trees, nests, or animal sightings, is crucial in understanding environmental patterns. This study compares MLGCP and SPIGPP, highlighting their strengths and weaknesses when prior knowledge of driving mechanisms is limited. Using synthetic and real datasets, we examine their predictive accuracy for the empirical K function and evaluate their performance under different interaction types and intensities.

The findings include:

MLGCP: Well-suited for capturing complex, unobserved spatial heterogeneities but struggles with strong repulsive interactions.
SPIGPP: Excels in estimating both the direction and magnitude of interactions, even under model misspecification.
This comparison offers valuable insights for researchers seeking appropriate models for analyzing multi-species point patterns in fields like ecology, epidemiology, and urban studies.

Installation

Dependencies
Ensure you have R (version 4.3.1 or above) installed on your system. Required R packages include those available in CRAN, along with the Multilogreg package (see below for instructions).



Pre-print can be found in;
Chathuri Samarasekara, Yan Wang, Ian Flint. A Comparison of Multivariate Log Gaussian Cox Process and Saturated Pairwise Interaction Gibbs Point Process. Authorea. September 13, 2024. DOI: 10.22541/au.172623817.72677555/v1

Please find that 
The original Multilogreg package can be accessed from https://github.com/kristianhessellund/Multilogreg.git and has compatibility and dependency issues. 

These have been fixed in https://github.com/IbTJensen/Multilogreg.git and can directly be installed from. 

The packages (Other than the ones available in CRAN) that needs for the codes to work can be installed using the code "Install required packages" in the repo.
