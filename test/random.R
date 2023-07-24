#!/usr/bin/env Rscript

# to use
# (echo sample,pheno; paste <(zcat ~/vcflib/samples/1kg-phaseIII-v5a.20130502.genotypes.chr22-16-16.5mb.vcf.gz | grep ^#CH | cut -f 10- | tr '\t' '\n') <(./random.R 2504 | tail -n+2 | awk '{ print $2 }') | tr '\t' ',') >x.csv

# Get the command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Get the number of random numbers to generate
n <- as.integer(args[1])

# Generate n random numbers between 0 and 1
random_numbers <- runif(n)

# Create a dataframe
df <- data.frame("RandomNumbers" = random_numbers)

# Print the dataframe
print(df)
