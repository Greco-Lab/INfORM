# INfORM
### A new tool for Inference of NetwOrk Response Modules

A novel computational method and its R and web-based implementations, to perform inference of gene network from transcriptome data and prioritization of key genes with central functional and topological role in the network.

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
