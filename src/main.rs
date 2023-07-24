use std::collections::HashMap;
use rust_htslib::bcf::{Reader, Read};
use csv::Reader as CsvReader;
use clap::Parser;

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

    // Write the header line
    print!("sample\tphenotype");
    for snp_id in snp_ids {
        print!("\t{}", snp_id);
    }
    println!();

    // Write each sample's genotype/phenotype data
    for result in phenotypes.records() {
        let record = result.unwrap();
        let sample = &record[0];
        let phenotype = &record[1];

        // convert sample string to &[u8]
        let sample_bytes = sample.as_bytes();
        if let Some(sample_genotypes) = sample_genotypes.get(sample_bytes) {
            print!("{}\t{}", sample, phenotype);
            for genotype in sample_genotypes {
                print!("\t{}", genotype);
            }
            println!();
        }
    }
}
