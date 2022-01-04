# INfORM
### A new tool for Inference of NetwOrk Response Modules

A novel computational method and its R and web-based implementations, to perform inference of gene network from transcriptome data and prioritization of key genes with central functional and topological role in the network.

Reference paper: Veer Singh Marwah, Pia Anneli Sofia Kinaret, Angela Serra, Giovanni Scala, Antti Lauerma, Vittorio Fortino, Dario Greco; INfORM: Inference of NetwOrk Response Modules, Bioinformatics, Volume 34, Issue 12, 15 June 2018, Pages 2136–2138

More information at: https://academic.oup.com/bioinformatics/article/34/12/2136/4841712

## Running the INfORM Docker image (suggested)

If you don't have docker installed on your system you can install it by following the instructions at  https://www.docker.com/get-docker.

The INfORM docker image is available at https://hub.docker.com/r/grecolab/inform


## Using INfORM source from GitHub


#### Install OS Dependencies
Pandoc - [https://pandoc.org/installing.html](https://pandoc.org/installing.html)  
Latex Compiling - [http://www.tug.org/texlive/](http://www.tug.org/texlive/)  
Perl - [https://www.perl.org/get.html](https://www.perl.org/get.html)  
Cairo - [https://cairographics.org/download/](https://cairographics.org/download/)  

#### Install R Dependencies
```R

  #Universal Bioconductor package installation function
  install.bioc <- function(pkg){
    vers <- getRversion()
    if (vers >= "3.6"){
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg)
    }else{
      if (!requireNamespace("BiocInstaller", quietly = TRUE)){
        source("https://bioconductor.org/biocLite.R")
        biocLite(pkg, suppressUpdates=TRUE)
      }else{
        BiocInstaller::biocLite(pkg, suppressUpdates=TRUE)
      }
    }
  }

  #Install CRAN dependencies
  cran_pkgs <- c("V8", "RSQLite", "TopKLists", "doParallel", "foreach", "igraph", "plyr", "shiny", "shinyjs", "shinyBS", "shinydashboard", "colourpicker", "DT", "R.utils", "treemap", "visNetwork", "abind", "radarchart", "randomcoloR", "Rserve", "WriteXLS", "gplots", "ggplot2", "devtools", "flock")
  cran_pkgs.inst <- cran_pkgs[!(cran_pkgs %in% rownames(installed.packages()))]
  if(length(cran_pkgs.inst)>0){
    print(paste0("Missing ", length(cran_pkgs.inst), " CRAN Packages:"))
    for(pkg in cran_pkgs.inst){
      print(paste0("Installing Package:'", pkg, "'..."))
      install.packages(pkg, repo="http://cran.rstudio.org", dependencies=TRUE)
      print("Installed!!!")
    }
  }

  #Install Bioconductor dependencies
  bioc_pkgs <- c("org.Hs.eg.db", "org.Mm.eg.db","org.Rn.eg.db", "GO.db", "AnnotationDbi", "GSEABase", "minet")
  bioc_pkgs.inst <- bioc_pkgs[!(bioc_pkgs %in% rownames(installed.packages()))]
  if(length(bioc_pkgs.inst)>0){
    print(paste0("Missing ", length(bioc_pkgs.inst), " Bioconductor Packages:"))
    for(pkg in bioc_pkgs.inst){
      print(paste0("Installing Package:'", pkg, "'..."))
      install.bioc(pkg)
      print("Installed!!!")
    }
  }

  #Install latest version of GOSemSim from GitHub
  print("Installing GOSemSim from GitHub!")
  devtools::install_github("GuangchuangYu/GOSemSim")
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

  # Start R session and run by using runApp()
  setwd("./INfORM_clone")
  library(shiny)
  runApp("INfORM-app/")
```
#### Dependencies and License
Please refer to the 'DESCRIPTION' and 'NAMESPACE' files for information about the license and dependencies required to run INfORM.
