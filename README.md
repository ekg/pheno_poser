# pheno_poser

This program takes a VCF file of genotypes and a CSV file of phenotypes as input. It performs the following steps:

1. Reads in the VCF and collects the sample names and genotype values for each variant, storing them in a matrix. 

2. Removes any genotype columns that are constant across all samples.

3. Transforms the genotype matrix so that rows are samples and columns are variants. 

4. Performs PCA on the genotype matrix, either using randomized PCA or full PCA depending on the `--full_pca` flag.

5. Takes the top N principal components as specified by `--pca_dims`.

6. Reads in the phenotype CSV file and extracts the phenotype values.

7. Builds a multiple linear regression model using the PCs as predictors and the phenotype as the response. 

8. Calculates the residuals for each sample from the fitted regression model.

9. Writes out a new file with the following columns:
   - Sample ID
   - Original phenotype value
   - Residual phenotype after regressing out PCs
   - PC1-PCN values
   - Genotype values

The residual phenotype represents the remaining phenotypic variance after accounting for population structure captured by the genotype PCs. This allows downstream analyses to find variants associated with phenotype independently of overall population structure.

The program takes the following arguments:

- `--vcf`: Input VCF file containing genotypes
- `--csv`: Input CSV file containing sample IDs and phenotype
- `--pca_dims`: Number of principal components to use in regression model  
- `--random_seed`: Random seed for randomized PCA
- `--oversampling`: Oversampling factor for randomized PCA
- `--full_pca`: Use full PCA instead of randomized PCA
