#include <iostream>
#include <stdexcept>
#include <string>
#include "htslib/vcf.h"
#include "libs/eigen-3.4.0/Dense"
#include "vcf.hpp"
using Eigen::ArrayXd;
using namespace std;


int count_sites_in_vcf(const string& vcf_path)
{
    // Prepare pointers
    htsFile* fp;
    bcf_hdr_t* hdr;
    bcf1_t* rec = bcf_init();

    // Initialise
    fp = vcf_open(vcf_path.c_str(), "r");
    hdr = vcf_hdr_read(fp);
    
    // Don't parse any sample information, for speed
    int r = bcf_hdr_set_samples(hdr, NULL, 0);
    if (r != 0) {
        throw std::invalid_argument("Failed to filter samples.");
    }

    // Count
    int n_sites = 0;
    while (vcf_read(fp, hdr, rec) == 0) {
        ++n_sites;
    }

    // Deallocate
    bcf_hdr_destroy(hdr);
    bcf_destroy(rec);
    vcf_close(fp);

    return n_sites;
}


ArrayXd calc_plafs_from_vcf(const string& vcf_path, double epsilon)
{
    // Prepare pointers
    htsFile* fp;
    bcf_hdr_t* hdr;
    bcf1_t* rec = bcf_init();

    // Initialise
    fp = vcf_open(vcf_path.c_str(), "r");
    hdr = vcf_hdr_read(fp);

    // Prepare arrays
    int n_samples = bcf_hdr_nsamples(hdr);
    int n_sites = count_sites_in_vcf(vcf_path);
    ArrayXd total_depth = ArrayXd::Constant(n_sites, 0);
    ArrayXd total_alts = ArrayXd::Constant(n_sites, 0);

    // Iterate over records
    int ix = 0;
    while(vcf_read(fp, hdr, rec) == 0) {
        // Store missingness for this record
        int n_samples_missing_data = 0;
        
        // Extract allelic depths
        int sz = 0;
        int* vals = NULL;
        if (bcf_get_format_int32(hdr, rec, "AD", &vals, &sz) > 0) {

            int N = 2 * n_samples;
            if (sz != N) {
                throw std::invalid_argument("Error in estimating PLAF from allelic depths. All sites must be bi-allelic.");
            }

            for (int j = 0; j < N; ++j) {
                
                // Check for missing data
                if (vals[j] == bcf_int32_missing) {
                    n_samples_missing_data += 1;
                    continue;
                }
                if (vals[j] == bcf_int32_vector_end) {
                    continue;
                }

                // Accumulate
                total_depth(ix) += vals[j];
                if (j % 2 == 1) {
                    total_alts(ix) += vals[j];
                }
            }

            if (n_samples - n_samples_missing_data < 10) {
                cout << "WARNING: Less than 10 samples have data for site " << ix << " in VCF." << endl;
            }

        } else {
            throw std::invalid_argument("Failed to load `AD` information. Exiting.");
        }

        ++ix;
        free(vals);
    }
    
    // Release memory
    bcf_hdr_destroy(hdr);
    bcf_destroy(rec);
    vcf_close(fp);

    ArrayXd plafs = total_alts / (total_depth + epsilon);
    return plafs;
}