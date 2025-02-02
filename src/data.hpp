#pragma once
#include <string>
#include <vector>
#include "htslib/vcf.h"
#include "libs/eigen-3.4.0/Dense"
#include "libs/rapidcsv.h"
#include "typedefs.hpp"
using namespace std;



/*
* Create `VCFData` structs from an input VCF file at `vcf_path`
* UNIMPLEMENTED:
* - Best to create this way to separate creation from use
* - However, want to be precise and efficient with memory management
* - May also want to allow for some polymorphism here
* NB: 
* - Do not copy this information
* - Could be big, allocate onto the heap
* - Probably best is with a shared pointer
*/
// class VCFDataFactory;



/*
* Encapsulate VCF data from `vcf_path` for a target sample, 
* given by `sample_name`
*
* TODO:
* - Implement some sort of factory for this
* - We really want the arrays to be const
* - We want this allocated on the heap (could be big)
* - We need to improve storage of chromosome information
*/
class VCFData
{
private:
    const string& vcf_path;             // Path to VCF
    const string& sample_name;          // Target sample name

    int calc_genome_length();
    
public:
    int n_samples;
    int n_sites;
    int n_sites_missing;
    int genome_length;                  // Length of genome in basepairs (bp)

    vector<string> chrom_names;
    ArrayXi chroms;
    ArrayXi pos;
    ArrayXi refs;                       // Read counts of REF allele
    ArrayXi alts;                       // Read counts of ALT allele
    ArrayXd plafs;                      // Population-level allele frequencies

    VCFData(const string& vcf_path, const string& sample_name);

    void check_valid() const;
    void print() const;
};



// Sure, I could use a template, but I only need these two types setups.
ArrayXd load_csv_as_double_array(const rapidcsv::Document& csv, const std::string& column_name);
ArrayXi load_csv_as_int_array(const rapidcsv::Document& csv, const std::string& column_name);


class InferredIBDPathData
{
private:
    const string& csv_path;             // Path to CSV file
    rapidcsv::Document csv;
public:
    std::vector<string> chrom_names;
    ArrayXi pos;
    ArrayXi refs;                       // Read counts of REF allele
    ArrayXi alts;                       // Read counts of ALT allele
    ArrayXd plafs;                      // Population-level allele frequencies
    ArrayXi ibd_path;                   // IBD state path
    int n_sites;

    InferredIBDPathData(const string& csv_path);
    void print() const;
};