#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include "vcf.hpp"
#include "data.hpp"
#include "constants.hpp"
using namespace std;


VCFData::VCFData(const string& vcf_path, const string& sample_name)
    : vcf_path(vcf_path),
    sample_name(sample_name),
    n_sites(count_sites_in_vcf(vcf_path)),
    n_sites_missing(0),
    chrom_names(n_sites),
    chroms(ArrayXi::Constant(n_sites, DEFAULT_INT)),
    pos(ArrayXi::Constant(n_sites, DEFAULT_INT)),
    refs(ArrayXi::Constant(n_sites, DEFAULT_INT)),
    alts(ArrayXi::Constant(n_sites, DEFAULT_INT)),
    plafs(ArrayXd::Constant(n_sites, DEFAULT_FLOAT))
{
    // Prepare pointers to file, header and record
    htsFile* fp;
    bcf_hdr_t* hdr;
    bcf1_t* rec = bcf_init();

    // Initialise
    fp = vcf_open(vcf_path.c_str(), "r");
    hdr = vcf_hdr_read(fp);

    // Count number of samples total
    n_samples = bcf_hdr_nsamples(hdr);

    // Filter to only target sample
    if (bcf_hdr_set_samples(hdr, sample_name.c_str(), 0) != 0) {
        throw std::invalid_argument("Error: could not filter to target sample.");
    }

    // Iterate over record and populate arrays
    int site_ix = 0;
    while (vcf_read(fp, hdr, rec) == 0) {
        
        // Site information
        chrom_names[site_ix] = bcf_hdr_id2name(hdr, rec->rid);
        chroms(site_ix) = rec->rid;
        pos(site_ix) = rec->pos;

        // Allelic depth information
        int sz = 0;
        int* vals = NULL;
        if (bcf_get_format_int32(hdr, rec, "AD", &vals, &sz) > 0) {
            
            if (sz == 1) { // could also be monoallelic, in which case treat as missing
                n_sites_missing += 1;
                refs(site_ix) = MISSING_AD_VALUE;
                alts(site_ix) = MISSING_AD_VALUE;
            } else if (sz != 2) {
                throw std::invalid_argument("Error parsing VCF: all sites must be bi-allelic.");
            } else {
                refs(site_ix) = vals[0];
                alts(site_ix) = vals[1];
            }

        } else {
            throw std::invalid_argument("Error: failed to read AD for a site.");
        }

        free(vals);
        ++site_ix;
    }

    // Release memory
    bcf_hdr_destroy(hdr);
    bcf_destroy(rec);
    vcf_close(fp);

    plafs = calc_plafs_from_vcf(vcf_path);
    genome_length = calc_genome_length();

    // Validity check
    check_valid();
}

void VCFData::check_valid() const
{
    if ((chroms == DEFAULT_INT).any()) {
        std::invalid_argument("Error parsing VCF!");
    }
    if ((pos == DEFAULT_INT).any()) {
        std::invalid_argument("Error parsing VCF!");
    }
    if ((refs == DEFAULT_INT).any()) {
        std::invalid_argument("Error parsing VCF!");
    }
    if ((alts == DEFAULT_INT).any()) {
        std::invalid_argument("Error parsing VCF!");
    }
    if ((plafs == DEFAULT_FLOAT).any()) {
        std::invalid_argument("Error parsing VCF!");
    }
}

int VCFData::calc_genome_length()
{

    // Initialise
    int i = 0;
    string chrom = chrom_names[i];
    int chrom_start = pos(i);
    int total_length = 0;
    ++i;

    // Iterate
    for (; i < n_sites; ++i) {

        // Continue if same chromosome
        if (chrom_names[i] == chrom) {
            continue;
        }

        // Else add length and update memory
        total_length += (pos(i-1) - chrom_start);
        chrom_start = pos(i);
        chrom = chrom_names[i];
    }

    // Terminate
    total_length += (pos(i-1) - chrom_start);

    return total_length;
}


void VCFData::print() const
{
    cout << "Loaded VCF data" << endl;
    cout << "  VCF path: " << vcf_path << endl;
    cout << "    No. sites: " << n_sites << endl;
    cout << "    No. samples: " << n_samples << endl;
    cout << "  Target sample: " << sample_name << endl;
    cout << "    No. sites with missing data for target sample: " << n_sites_missing << endl;

    // TODO: should make a macro or inline function for these...
    cout << "  Site data:" << endl;
    cout << "    Chromosomes: " << chrom_names[0] << " ... " << chrom_names[n_sites - 1] << endl;
    cout << "    Positions: " << pos(0) << " ... " << pos(n_sites - 1) << endl;
    cout << "    REFs: " << refs(0) << " ... " << refs(n_sites - 1) << endl;
    cout << "    ALTs: " << alts(0) << " ... " << alts(n_sites - 1) << endl;
    cout << "    PLAFs: " << plafs(0) << " ... " << plafs(n_sites - 1) << endl;
}

