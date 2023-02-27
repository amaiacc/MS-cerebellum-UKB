# following: https://github.com/josefin-werme/LAVA

## prerequisites
install.packages("BiocManager")
BiocManager::install("snpStats")
install.packages("devtools")

# install lava
devtools::install_github("https://github.com/josefin-werme/LAVA.git")

# test
library(LAVA)
