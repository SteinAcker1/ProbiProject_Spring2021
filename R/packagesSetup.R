# Installing from Bioconductor
if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos="http://cran.us.r-project.org", quiet = TRUE)
}
BiocManager::install(version = "3.13")
BiocManager::install("dada2")

# Installing from CRAN
if(!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos="http://cran.us.r-project.org", quiet = TRUE)
}
devtools::install_version("tidyverse", version = "1.3.1", repos="http://cran.us.r-project.org", quiet = TRUE)
devtools::install_version("vegan", version = "2.5.7", repos="http://cran.us.r-project.org", quiet = TRUE)
devtools::install_version("ggrepel", version = "0.9.1", repos="http://cran.us.r-project.org", quiet = TRUE)
devtools::install_version("viridis", version = "0.6.1", repos="http://cran.us.r-project.org", quiet = TRUE)
devtools::install_version("tidyMicro", version = "1.47", repos="http://cran.us.r-project.org", quiet = TRUE)
devtools::install_version("latex2exp", version = "0.5.0", repos="http://cran.us.r-project.org", quiet = TRUE)

