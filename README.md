# INfORM
### A new tool for Inference of NetwOrk Response Modules

A novel computational method and its R and web-based implementations, to perform inference of gene network from transcriptome data and prioritization of key genes with central functional and topological role in the network.

#### Install Dependencies
```R
  #Install CRAN dependencies
  cran_pkgs <- c("V8", "RSQLite", "TopKLists", "doParallel", "foreach", "igraph", "plyr", "shiny", "shinyjs", "shinyBS", "shinydashboard", "colourpicker", "DT", "R.utils", "treemap", "visNetwork", "abind")
  install.packages(cran_pkgs, repo="http://cran.rstudio.org", dependencies=T)
  
  #Install Bioconductor dependencies
  source("http://bioconductor.org/biocLite.R")
  bioc_pkgs <- c("org.Hs.eg.db", "org.Mm.eg.db", "GO.db", "AnnotationDbi", "GSEABase", "minet", "GOSemSim")
  biocLite(bioc_pkg, suppressUpdates=T)
```

#### How to run INfORM from GitHub
```R
  # Load 'shiny' library
  library(shiny)

  # Using runGitHub
  runGitHub("INfORM", "Greco-Lab", subdir="INfORM-app")

  # Using the archived file
  runUrl("https://github.com/Greco-Lab/INfORM/archive/master.tar.gz", subdir="INfORM-app")
  runUrl("https://github.com/Greco-Lab/INfORM/archive/master.zip", subdir="INfORM-app")
```

#### How to run locally
```R
  # Clone the git repository
  git clone https://github.com/Greco-Lab/INfORM INfORM_clone

  # Run by using runApp()
  setwd("~/INfORM_clone")
  runApp("INfORM-app/")
```
#### Dependencies and License
Please refer to the 'DESCRIPTION' and 'NAMESPACE' files for information about the license and dependencies required to run INfORM.
