use std::collections::{HashMap, HashSet};
use rust_htslib::bcf::{Reader, Read, record};
use csv::Reader as CsvReader;
use clap::Parser;
//use statrs::statistics::PCA;
use pca::PCA;
use ndarray::{Array2, s};
//use ndarray::arr2;
//use petal_decomposition::{Pca, PcaBuilder, RandomizedPcaBuilder};
use linregress::{FormulaRegressionBuilder, RegressionDataBuilder};


#[derive(Parser)]
#[command(name = "Genotype Phenotype Mapper")]
#[command(author = "Your Name <your.email@gmail.com>")]
#[command(version = "0.1")]
#[command(about = "Maps genotype to phenotype based on input VCF and phenotype CSV files.", long_about = None)]
struct Cli {
    /// Input VCF file
    #[arg(short, long)]
    vcf_in: String,
    /// Input phenotype CSV file
    #[arg(short, long)]
    csv_in: String,
    /// Keep this many dimensions of PCA
    #[arg(short, long, default_value = "10")]
    pca_dims: usize,
}

fn main() {
    let cli = Cli::parse();

    // Open the VCF file
    let mut vcf = Reader::from_path(&cli.vcf_in).unwrap();
    // collect sample names
    let header = vcf.header().clone();
    let sample_names = header.samples().to_vec();
    //let sample_names = sample_binding.iter().map(|s| s.to_owned());
    
    let mut sample_genotypes: HashMap<&[u8], Vec<i32>> = HashMap::new();
    for sample in sample_names.clone() {
        sample_genotypes.insert(sample, Vec::new());
    }

    let mut snp_ids: Vec<String> = Vec::new();

    for (_i, record_result) in vcf.records().enumerate() {
        let record = record_result.expect("Fail to read record");
        let mut s = String::new();
        // get SNP Id
        snp_ids.push(String::from_utf8(record.id().to_vec()).unwrap());
        for allele in record.alleles() {
            for c in allele {
                s.push(char::from(*c))
            }
            s.push(' ')
        }
        // 0-based position and the list of alleles
        //println!("Locus: {}, Alleles: {}", record.pos(), s);
        // number of sample in the vcf
        let sample_count = usize::try_from(record.sample_count()).unwrap();

        // Counting ref, alt and missing alleles for each sample
        let gts = record.genotypes().expect("Error reading genotypes");
        for sample_index in 0..sample_count {
            // for each sample
            //let mut n_ref = vec![0; sample_count];
            //let mut n_alt = vec![0; sample_count];
            //let mut n_missing = vec![0; sample_count];
            let mut _n_ref = 0;
            let mut n_alt = 0;
            let mut _n_missing = 0;
            for gta in gts.get(sample_index).iter() {
                // for each allele
                match gta.index() {
                    Some(0) => _n_ref += 1,  // reference allele
                    Some(_) => n_alt += 1,  // alt allele
                    None => _n_missing += 1, // missing allele
                }
            }
            let s = sample_names.get(sample_index).unwrap();
            sample_genotypes.get_mut(s).unwrap().push(n_alt);
        }
    }


    // Open the phenotype CSV file
    let mut phenotypes = CsvReader::from_path(&cli.csv_in).unwrap();


    // Remove constant genotype columns
    eprintln!("Removing constant genotype columns...");
    // for each column in the matrix
    // for each sample, get the column value
    // track if they are all the same using boolean logic

    // first, we get the length of the genotype vectors from the first sample
    let mut sample_genotype_lengths = Vec::new();
    for (_sample, genotype) in sample_genotypes.iter() {
        sample_genotype_lengths.push(genotype.len());
    }
    // assert they are all equal in length
    assert!(sample_genotype_lengths.iter().all(|&item| item == sample_genotype_lengths[0]));
    // take the first one's length
    let genotypes_length = sample_genotype_lengths[0];
    eprintln!("Genotype length: {}", genotypes_length);

    // now from 0 to genotypes_length, we look across samples and check if they are all the same
    // we will record the constant columns here:
    let mut constant_columns = HashSet::new();
    for i in 0..genotypes_length {
        let mut is_constant = true;
        let mut last_value = 0;
        // iterate over all samples, accessing the genotype at this offset
        for j in 0..sample_names.len() {
            let sample = sample_names.get(j).unwrap();
            let genotype = sample_genotypes.get(sample).unwrap();
            let value = genotype.get(i).unwrap();
            if j == 0 {
                last_value = *value;
            } else if last_value != *value {
                is_constant = false;
                break;
            }
        }
        if is_constant {
            constant_columns.insert(i);
        }
    }

    let mut filtered_genotypes = HashMap::new();
    for (sample, genotype) in sample_genotypes.iter() {
        let mut filtered_genotype = Vec::new();
        for (i, g) in genotype.iter().enumerate() {
            if !constant_columns.contains(&i) {
                filtered_genotype.push(*g);
            }
        }
        filtered_genotypes.insert(sample, filtered_genotype);
    }

    eprintln!("Number of samples: {}", filtered_genotypes.len());
    eprintln!("Number of SNPs: {}", filtered_genotypes.values().next().unwrap().len());

    let mut pca_data = Vec::new();
    for (_sample, genotype) in sample_genotypes.iter() {
        let mut pca_sample = Vec::new();
        for (i, g) in genotype.iter().enumerate() {
            if !constant_columns.contains(&i) {
                pca_sample.push(*g as f64);
            }
        }
        pca_data.push(pca_sample);
    }

    // Check for NaN or infinite values
    for (i, row) in pca_data.iter().enumerate() {
        if row.iter().any(|&x| x.is_nan() || x.is_infinite()) {
            eprintln!("Row {} contains NaN or infinite values!", i);
            return; // or handle the error as you prefer
        }
    }

    // save the array to disk in a TSV file
    use std::fs::File;
    use std::io::Write;
    let mut pca_data_file = File::create("pca_data.tsv").unwrap();
    for row in pca_data.iter() {
        for (i, value) in row.iter().enumerate() {
            if i > 0 {
                pca_data_file.write_all(b"\t").unwrap();
            }
            pca_data_file.write_all(format!("{}", value).as_bytes()).unwrap();
        }
        pca_data_file.write_all(b"\n").unwrap();
    }

    eprintln!("Building input array...");
    // build the array, remembering that we have things in the format of rows of samples (each a vector) and columns of SNPs
    let x = Array2::from_shape_vec((pca_data.len(), pca_data[0].len()), pca_data.into_iter().flatten().collect()).unwrap();

    eprintln!("First few rows of x:");
    for row in x.slice(s![..3, ..]).rows() {
        eprintln!("{:?}", row);
    }
    
    //let x = Array2::from_shape_vec((pca_data[0].len(), pca_data.len()), pca_data.iter().flatten().cloned().collect()).unwrap();

    eprintln!("Performing PCA...");
    let mut pca = PCA::new();
    pca.fit(x.clone(), None).unwrap();

    let components = pca.transform(x).unwrap();
    // print dimensions
    eprintln!("Components shape: {:?}", components.shape());
    // slice to the first cli.pca_dims
    let components = components.slice(s![.., ..cli.pca_dims]).to_owned();
    eprintln!("Components sliced shape: {:?}", components.shape());

    // print a slice of the components for debugging
    //eprintln!("{:?}", components.slice(s![..10, ..]));
    
    //let mut pca = RandomizedPcaBuilder::new(cli.pca_dims).build(); // Keep N dimensions.
    //let mut pca = PcaBuilder::new(cli.pca_dims).build(); // Keep N dimensions.
    //pca.fit(&x).unwrap();

    //eprintln!("Explained variance: {:?}", pca.explained_variance_ratio());
    //let _s = pca.singular_values();            // [2_f64, 0_f64]
    //let _v = pca.explained_variance_ratio();   // [1_f64, 0_f64]
    //let _y = pca.transform(&x).unwrap();       // [-2_f64.sqrt(), 0_f64, 2_f64.sqrt()]

    eprintln!("Collecting components...");
    //let components = pca.components();

    // Build GLM model
    // follow this model let data = vec![("Y", y), ("X1", x1), ("X2", x2), ("X3", x3)];
    // we're going to build a model with the PCs in components as X and the phenotypes as Y
    let mut data = Vec::new();
    // add the phenotype
    let mut y = Vec::new();
    // collect phenotype records
    let mut phenotype_records = Vec::new();
    for record in phenotypes.records() {
        phenotype_records.push(record.unwrap());
    }

    // Write each sample's genotype/phenotype data
    for record in phenotype_records.clone().into_iter() {
        let phenotype = record.get(1).unwrap().parse::<f64>().unwrap();
        y.push(phenotype);
    }
    data.push((format!("Y"), y));
    for (i, component) in components.rows().into_iter().enumerate() {
        let mut component_data = Vec::new();
        for (j, value) in component.iter().enumerate() {
            component_data.push(*value);
        }
        eprintln!("Component {}: {:?}", i, component_data.len());
        //eprintln!("Component {}: {:?}", i, component_data);
        data.push((format!("PC{}", i+1), component_data));
    }

    let pcs = (1..cli.pca_dims+1).map(|i| format!("PC{}", i)).collect::<Vec<String>>();

    for (i, pc) in pcs.iter().enumerate() {
        eprintln!("{}, {}, {}", i, pc, data[i+1].1.len());
    }

    // Write the header line
    print!("sample\tphenotype\tresidual");
    // Write columns for each PC
    for i in 0..cli.pca_dims {
        print!("\tPC{}", i+1);
    }
    for snp_id in snp_ids {
        print!("\t{}", snp_id);
    }
    println!();
    
    let mut sample_phenotype = HashMap::new();
    for record in phenotype_records {
        let sample = String::from(record.get(0).unwrap());
        let phenotype = String::from(record.get(1).unwrap());
        sample_phenotype.insert(sample, phenotype);
    }

    for (i, (sample_bytes, genotypes)) in sample_genotypes.iter().enumerate() {
        // convert &[u8] to string
        let sample = String::from_utf8_lossy(sample_bytes).into_owned();
        if let Some(phenotype) = sample_phenotype.get(&sample) {
            eprintln!("Sample {} has phenotype {}", sample, phenotype);
            print!("{}\t{}\t{}", sample, phenotype, 0);
            // print out the PCs
            for j in 0..cli.pca_dims {
                print!("\t{}", components.row(i)[j]);
            }
            for genotype in genotypes {
                print!("\t{}", genotype);
            }
            println!();
        }
    }

    eprintln!("Done!");

    /*
    eprintln!("Building GLM model...");

    let data = RegressionDataBuilder::new().build_from(data).unwrap();        
    //let formula = "Y ~ X1 + X2 + X3";
    let model = FormulaRegressionBuilder::new().data(&data).data_columns("Y", pcs).fit().unwrap();
    let _parameters: Vec<_> = model.iter_parameter_pairs().collect();
    let _pvalues: Vec<_> = model.iter_p_value_pairs().collect();
    let _standard_errors: Vec<_> = model.iter_se_pairs().collect();

    // get the residuals for each sample
    let residuals = model.residuals();
    eprintln!("Residuals: {:?}", residuals.len());
    
    // Write the header line
    print!("sample\tphenotype\tresidual");
    // Write columns for each PC
    for i in 0..cli.pca_dims {
        print!("\tPC{}", i+1);
    }
    for snp_id in snp_ids {
        print!("\t{}", snp_id);
    }
    println!();

    //println!("sample genotypes size {}", sample_genotypes.len());
    // Write each sample's genotype/phenotype data
    let mut sample_phenotype = HashMap::new();
    for record in phenotype_records {
        let sample = String::from(record.get(0).unwrap());
        let phenotype = String::from(record.get(1).unwrap());
        sample_phenotype.insert(sample, phenotype);
    }

    for (i, (sample_bytes, genotypes)) in sample_genotypes.iter().enumerate() {
        // convert &[u8] to string
        let sample = String::from_utf8_lossy(sample_bytes).into_owned();
        if let Some(phenotype) = sample_phenotype.get(&sample) {
            print!("{}\t{}\t{}", sample, phenotype, residuals[i]);
            // print out the PCs
            for j in 0..cli.pca_dims {
                print!("\t{}", components.row(j)[i]);
            }
            for genotype in genotypes {
                print!("\t{}", genotype);
            }
            println!();
        }
}
    */

}
