# ~~~ _init.R ~~~
# * ONLY RUN THIS SCRIPT TO GENERATE A NEW RENV (not recommended)
# - 1 ~ renv (virtual environment)
setwd('D:/repos/baylor/tlgll-seurat-characterization')

install.packages('renv')
library(renv)
renv::init()

# - 2 ~ CRAN

# - 3 ~ Bioconductor

# - 4 ~ Snapshots
renv::snapshot()
writeLines(capture.output(sessionInfo()), 'session_info.txt')