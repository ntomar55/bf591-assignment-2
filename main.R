#!/usr/bin/Rscript
## Author: Taylor Falk
## tfalk@bu.edu
## BU BF591
## Assignment Week 2

#### Bioconductor ####
load_bioconductor <- function() {
  # have students find correct install lines from https://bioconductor.org/install/ ?
  # note that this is different for the current version of Bioconductor than the 528 tutorial
  if (!require("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install(c('affy','affyPLM','sva',
                         'AnnotationDbi','hgu133plus2.db'), 
                       version = "3.13")
  library(affy)
  library(affyPLM)
  library(sva)
  library(AnnotationDbi)
  library(hgu133plus2.db)
}

load_expression <- function(filepath) {
  # Load the csv of expression data from the filepath variable, returns stored object
  expression <- read.csv(filepath, header = TRUE, sep = " ")
  return(expression)
}

#### Data types ####


#### Markdown ####
#source(knitr::purl('report.Rmd'))

if(interactive()) {
  # load_bioconductor()
  # expr <- load_expression('/project/bf528/project_1/data/example_intensity_data.csv')
}
