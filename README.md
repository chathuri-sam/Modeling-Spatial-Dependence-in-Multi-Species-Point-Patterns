# A Comparison of Multivariate Log Gaussian Cox Process and Saturated Pairwise Interaction Gibbs Point Process

The study of the spatial point patterns in ecology, such as the records of the observed locations of trees, shrubs, nests, burrows, or documented animal presence, relies on multivariate point process models. This study aims to compare the efficacy and applicability of two prominent multivariate point process models, the multivariate log Gaussian Cox process (MLGCP) and the Saturated Pairwise Interaction Gibbs Point Process model (SPIGPP) , highlighting their respective strengths and weaknesses in various scenarios. Using synthetic and real datasets, we assessed both models based on their predictive accuracy of the empirical K function. Our analysis revealed that both MLGCP and SPIGPP effectively identify and capture mild to moderate attractions and regulations. MLGCP struggles to capture repulsive associations as they intensify. In contrast, SPIGPP can well estimates both the direction and magnitude of interactions even when the model is miss-specified. Both models present unique advantages: MLGCP is particularly effective when there is a need to account for complex, unobserved heterogeneities that vary across space, while SPIGPP is suitable when interactions between points are the primary focus. The choice between these models should be guided by the specific needs of the research question and data characteristics.

Pre-print can be found in;
Chathuri Samarasekara, Yan Wang, Ian Flint. A Comparison of Multivariate Log Gaussian Cox Process and Saturated Pairwise Interaction Gibbs Point Process. Authorea. September 13, 2024. DOI: 10.22541/au.172623817.72677555/v1


The original Multilogreg package can be accessed from https://github.com/kristianhessellund/Multilogreg.git and has compatibility and dependency issues. 

These have been fixed in https://github.com/IbTJensen/Multilogreg.git and can directly be installed from. 
