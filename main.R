#!/usr/bin/Rscript
## Author: Taylor Falk
## tfalk@bu.edu
## BU BF591
## Assignment Week 2

#### Bioconductor ####
# it is standard among R packages to define libraries and packages at the 
# beginning of a script. Also note that a package should NOT be installed every 
# time a script runs.
# The bioconductor repository has installation instructions for biomaRt: 
# https://bioconductor.org/install/

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
if (!require("biomaRt", quietly = TRUE)){
  BiocManager::install("biomaRt")
}
library(biomaRt)
library(tidyverse)

load_expression <- function(filepath) {
  # Load the csv of expression data from the filepath variable, returns a tibble
  return(read_csv(filepath))
}

#' Title
#'
#' @param tibble A tibble of expression values, rows by probe and columns by sample.
#'
#' @return A vector of affymetrix probe names from the input expression data tibble. 
#' These names match the rows with 15% or more of the expression values about log2(15)
#' @export
#'
#' @examples
filter_15 <- function(tibble){
  # for each gene, at least 15% of the gene-expression values must be > log2(15)
  # return a list of sample names with 15% of values greater than log2(15)
  percent_gt <- function(row) {
    # functions can be defined inside other functions, which can be one style to
    # make your code more repeatable
    boolean_row <- row > log2(15)
    percent_row <- sum(boolean_row)/length(row)
    return(percent_row)
  }
  # there are many ways in R to apply our function to the entire tibble
  # fastest would be an lapply, but a for loop would work (just slowly)
  row_pct <- apply(tibble[2:ncol(tibble)], percent_gt, MARGIN = 1) # don't capture first row
  boolean_rows <- which(row_pct > 0.15)
  return(tibble[boolean_rows, 1])
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

#### ggplot ####

plot_ggplot <- function(dataframe) {
  # starting with a dataframe, plot a distribution of all of the expression
  # data points. return a ggplot object (i.e. if you run this function, it should
  # plot a graph in rstudio)
  long_df <- gather(expr)
  p <- ggplot(long_df, aes(x=value)) +
    geom_histogram() +
    geom_vline(xintercept = log2(15))
  return(p)
}

#### gene name conversion ####

#' Convert affymetrix array names into hgnc_symbol IDs using biomaRt
#'
#' @param affy_vector 
#'
#' @return 
#' @export
#'
#' @examples
affy_to_hgnc <- function(affy_vector) {
  ensembl <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  newNames <- getBM(attributes = c('affy_hg_u133_plus_2', 'hgnc_symbol'), 
        filters = c('affy_hg_u133_plus_2'),
        values = testMart[1:100], 
        mart = ensembl)
  return(newNames)
}

#### Markdown ####
#source(knitr::purl('report.Rmd'))

# this if statement means that this segment of code will only run when in 
# interactive mode, i.e. RStudio. You can include examples for your functions 
# here, and not worry about loading a lot of data when you source() this from
# another script or RMarkdown file. 
# This kind of setup may also prove useful if you select "source on save", which 
# runs your entire script every time it saves.
if(interactive()) {
  tmp_csv <- read.csv('/project/bf528/project_1/data/example_intensity_data.csv',
                      header = T, sep = " ")
  tmp_csv <- cbind(probeids = row.names(tmp_csv), tmp_csv)
  write_csv(tmp_csv, file = "temp_intensity_data.csv")
  expr <- load_expression("temp_intensity_data.csv")
  samples <- filter_15(expr)
  #p <- plot_ggplot(expr)
}
