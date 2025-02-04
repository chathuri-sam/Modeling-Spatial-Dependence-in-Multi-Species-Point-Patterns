# Clear environment
remove(list = ls())

# Install required packages from CRAN
install.packages(c("devtools", "spatstat", "ggplot2","ecespa")) 

# Install other dependencies, that are not in CRAN
library(devtools)

# Install RandomFieldsUtils which is not on CRAN anymore and fails to compile on some configs
# when installing an old version
devtools::install_github("iflint1/RandomFieldsUtils")

# Same as above with RandomFields package
devtools::install_github("iflint1/RandomFields")

# Old version of geostatsp that works with old Matrix package
devtools::install_version("geostatsp", dependencies = TRUE, version = "2.0.0", repos = "http://cran.us.r-project.org")

# Install Multilogreg package
devtools::install_github("chathuri-sam/Multilogreg_updated")

# Install ppjsdm package 
devtools::install_github("iflint1/ppjsdm")
