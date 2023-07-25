# vcf_pheno_poser

This program takes a VCF file of genotypes and a CSV file of phenotypes as input. It performs the following steps:

1. Filters the VCF to remove genotypes that are constant across all samples.

2. Transposes the genotype matrix so samples are columns and variants are rows.

3. Performs PCA on the genotype matrix and takes the top N principal components.

4. Fits a linear regression model using the principal components as predictors and the phenotype as the response. 

5. Calculates the residuals from the regression model for each sample.

6. Writes out a new file with the following columns:
   - Sample ID
   - Original phenotype
   - Residual PCA-adjusted phenotype
   - Genotypes

The purpose of this procedure is to incorporate population structure information into the phenotype by modeling it using the top principal components from the genotype matrix. The residuals represent the remaining phenotypic variance after accounting for population structure. This allows downstream analyses to find variants associated with phenotype independently of overall population structure.

The program takes the following arguments:

- `--vcf`: Input VCF file containing genotypes
- `--csv`: Input CSV file containing sample IDs and phenotype 
- `--pca_dims`: Number of principal components to use in regression model

The output file will contain the sample ID, original phenotype, residual phenotype, and genotype values for sites in the filtered VCF. This allows users to test for genotype-phenotype associations using the residual phenotype as the outcome.
