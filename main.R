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
  BiocManager::install(c('affy','affyPLM','sva', 'Biobase',
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

filter_15 <- function(dataframe){
  # for each gene, at least 15% of the gene-expression values must be > log2(15)
  # return a list of sample positions with 15% of values greater than log2(15)
  percent_gt <- function(row) {
    # functions can be defined inside other functions, which can be one style to
    # make your code more repeatable
    boolean_row <- row[1,] > log2(15)
    percent_row <- sum(boolean_row)/length(row)
    return(percent_row)
  }
  # there are many ways in R to apply our function to the entire data frame
  # fastest would be an apply, but we can use a for loop for now
  output_list <- data.frame(samples = NULL)
  for (row in seq(1, length(row.names(dataframe)))) {
    if (percent_gt(dataframe[row,]) >= 0.15) {
      output_list <- rbind(output_list, row)
    }
  }
  return(output_list)
}

#### Data types ####
join_str <- function(string1, string2) {
  # take parameters string1 and string2, both strings, and return one combined string
  # join_str("hello ", "world") returns  "hello world"
  # NOTE: could also do something more complicated, a string search maybe
  return(paste0(string1, string2))
}

boolean_filter <- function(boolean_array, int_array){
  # given a vector of booleans and an equally sized vector of integers, return 
  # a vector of all the integers corresponding to TRUE positions in the boolean.
  # e.g. boolean_filter(c(TRUE, FALSE, TRUE), c(1, 2, 3)) would return [1] 1 3
  return(int_array[boolean_array])
}

df_trim <- function(df, row_name, col_name) {
  # given a data.frame df, and a vector of row_names and col_names, return a new
  # data.frame with _only_ those rows and columns selected.
  return(df[row_name, col_name])
}

row.names(df) <- c("row1", "row2", "row3")
df_trim(df, c("row2", "row3"), c("a"))

#### Markdown ####
#source(knitr::purl('report.Rmd'))

if(interactive()) {
  # load_bioconductor()
  # expr <- load_expression('/project/bf528/project_1/data/example_intensity_data.csv')
  # samples <- filter_15(expr)
}
