# T-LGLL Seurat Characterization
This repository contains my final project for **BINF 3360**, focused on characterizing how Seurat parameter choices 
affect single-cell RNA-seq analysis results in T-cell large granular lymphocytic leukemia (T-LGLL).

## Reproducibility
This project uses `renv` for dependency management.

To reproduce the R environment:
```r
install.packages("renv")
renv::restore()
```

R version used: 4.5.3

## Workflow
This project has been organized in a sequential file workflow in the `scripts/` directory. 
It is strongly encouraged to run each file in its intended order, from 00_config to 05_figures.

## Dataset
- **GEO Accession ID:** [GSE168859](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168859)