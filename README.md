# virtimepredicteR
<!-- badges: start -->
<!-- badges: end -->

The goal of virtimepredicteR is to:
- Take an in-frame fasta file alignment (DNAStringSet object) and 
- Calculate Average Pairwise Distance (APD on all nucleotides)
- Also calculate APD on the 3rd codon position 
- Select best model for that gene region, and 
- Predict viremic time

For each gene region, model options include:
- Linear Regression of viremic time (y) against APD (x)
- Linear Regression of y against APD from 3rd codon positions
- Piecewise Linear Regression of y against on APD
- Piecewise Linear Regression of y against APD from 3rd codon positions

## Installation

Installing the development version of virtimepredicteR:

``` r
if (!require("devtools", character.only = TRUE)){install.packages("devtools")} else {require("devtools")}
devtools::install_github("ekankaka/virtimepredicteR")
```

## Example

This is a basic example which shows you how to solve a common problem:

## load the package
``` r
library(virtimepredicteR)
```
## load your in-frame fasta file alignment using the Biostrings package)
Here, we are using the example fasta alignment which comes with the package
``` r
fasta_path <- system.file("extdata", "example_GAG_P17.fasta", package = "virtimepredicteR")
myfasta <- Biostrings::readDNAStringSet(fasta_path)
```


## gene region
```r
# Should be one of the regions for which we have good predictive models. 
# These include: "GAG_P17", "POL_RT", "POL_IN", "GAG_P24", "POL_RT_TRIMMED", 
# "VPU", "REV1", "GP120_C1", "REV" 

this_region <- "GAG_P17"
```
## Calculate Average Pairwise Diversity 
APD on all nucleotides (APD)
```r
APD <- get_APD(fas = myfasta)[1]
```

APD on 3rd codon positions of in-frame alignment
```r
APD_codons <- get_APD(fas = myfasta)[2]
```

## APD version to use 
Determine this from the best_models dataset that comes with the package.
```r
APD_version_to_use = best_models$APD_version[best_models$gene_region = this_region]

APD_val = case_when(
  APD_version_to_use == "APD" ~ APD,
  APD_version_to_use == "APD_codons" ~ APD_codons)
  
APD_val
```

## Estimate viremic time
```r
years_untreated <- get_viremictime(APD_value = APD_val, gene_region = this_region)
```

## View result
```r
years_untreated
```
