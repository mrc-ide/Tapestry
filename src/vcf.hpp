#pragma once
#include <stdexcept>
#include <string>
#include "htslib/vcf.h"
#include "libs/eigen-3.4.0/Dense"
using Eigen::ArrayXd;
using namespace std;

/*
* Count the number of loci / variants in the VCF file at `vcf_path`
* Here, equivalent to number of rows
*/
int count_sites_in_vcf(const string& vcf_path);


/*
* Calculate population-level allele frequencies (PLAFs) for all sites
* in the VCF file at `vcf_path`
* NB: 
* - Computed from allelic depths
* - `epsilon` keeps estimates valid for sites with zero depth
* TODO:
* - Add checks for file opening
* - Add checks for AD field
*/
ArrayXd calc_plafs_from_vcf(const string& vcf_path, double epsilon = 0.001);
