#include <iostream>
#include <string>
#include <vector>
#include "htslib/vcf.h"
#include "libs/eigen-3.4.0/Dense"
using Eigen::ArrayXd;
using Eigen::ArrayXXi;
using namespace std;

// Thoughts:
// - Probably better to keep data in Eigen, for consistency
// - For chromosome, can have chromosome name and integer ID

struct Data
{
    int n_samples;
    int n_loci;
    vector<string> chroms;
    vector<int> pos;
    Data() {};
};

int calc_num_loci_from_vcf(const string& vcf_path)
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
    int n_loci = 0;
    while (vcf_read(fp, hdr, rec) == 0) {
        ++n_loci;
    }

    // Deallocate
    bcf_hdr_destroy(hdr);
    bcf_destroy(rec);
    vcf_close(fp);

    return n_loci;
}

ArrayXd calc_plafs_from_vcf(const string& vcf_path)
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
    int n_loci = calc_num_loci_from_vcf(vcf_path);
    ArrayXd total_depth(n_loci);
    ArrayXd total_alts(n_loci);

    // Iterate over records
    int ix = 0;
    while(vcf_read(fp, hdr, rec) == 0) {
        
        // Extract allelic depths
        int sz = 0;
        int* vals = NULL;
        if (bcf_get_format_int32(hdr, rec, "AD", &vals, &sz) > 0) {
            if (sz != 2*n_samples) {
                throw std::invalid_argument("Error. All sites must be bi-allelic.");
            }

            // Store total depth, and alternative depth
            for (int j = 0; j < n_samples; ++j) {
                total_depth(ix) += vals[j];
                if (j % 2 == 1) {
                    total_alts(ix) += vals[j];
                }
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

    ArrayXd plafs = total_alts / (total_depth + 0.001);
    return plafs;
}




int main(int argc, char* argv[])
{
    if (argc != 3) {
        cerr << "Usage: ./a.out <input_file.vcf> <sample_name>" << endl;
        exit(1);
    }
    string input_vcf = argv[1];
    string sample_name = argv[2];
    cout << "Input VCF: " << input_vcf << endl;
    cout << "Sample name: " << sample_name << endl;

    // Prepare pointers to file, header and record
    htsFile* fp;
    bcf_hdr_t* hdr;
    bcf1_t* rec = bcf_init();

    // Initialise
    fp = vcf_open(argv[1], "r");
    hdr = vcf_hdr_read(fp);

    // Here we filter to target sample
    bool filter_sample = true;
    if (filter_sample) {
        cout << "Filtering to target sample..." << endl;
        int r = bcf_hdr_set_samples(hdr, sample_name.c_str(), 0);
        if (r == 0) {
            cout << "Completed successfully." << endl;
        } else {
            cout << "Failed! Exit code: " << r << endl;
            exit(1);
        }
    }

    // Prepare storage
    int n_samples = bcf_hdr_nsamples(hdr);
    vector<string> chrom_names;
    vector<int> chroms; // these are a bit weird
    vector<int> pos;
    vector<int> refs;
    vector<int> alts;
    vector<double> wsafs;

    // Iterate
    while (vcf_read(fp, hdr, rec) == 0) {

        // Variant related information
        chrom_names.push_back(bcf_hdr_id2name(hdr, rec->rid));
        chroms.push_back(rec->rid);
        pos.push_back(rec->pos);

        // Allelic depths
        int sz = 0;
        int* vals = NULL;
        if (bcf_get_format_int32(hdr, rec, "AD", &vals, &sz) > 0) {
            if (sz != 2) {
                throw std::invalid_argument("Error. All sites must be bi-allelic.");
            }
            refs.push_back(vals[0]);
            alts.push_back(vals[1]);
            double wsaf = vals[1] / (vals[0] + vals[1] + 0.001);
            wsafs.push_back(wsaf);
        }

        // Release memory
        free(vals);
    }

    // Release memory
    bcf_hdr_destroy(hdr);
    bcf_destroy(rec);
    vcf_close(fp);

    // Compute number of loci
    int n_loci = chroms.size();

    // Write results
    cout << "Data" << endl;
    cout << "  No. Samples: " << n_samples << endl;
    cout << "  No. Loci: " << n_loci << endl;
    cout << "  Chrom. names: " << chrom_names[0] << " ... " << chrom_names[n_loci - 1] << endl;
    cout << "  Chrom. IDs: " << chroms[0] << " ... " << chroms[n_loci - 1] << endl;
    cout << "  Positions: " << pos[0] << " ... " << pos[n_loci - 1] << endl;
    cout << "  REFs: " << refs[0] << " ... " << refs[n_loci - 1] << endl;
    cout << "  ALTs: " << alts[0] << " ... " << alts[n_loci - 1] << endl;
    cout << "  WSAFs: " << wsafs[0] << " ... " << wsafs[n_loci - 1] << endl;
    cout << "Done." << endl;

    // Compute PLAFS
    cout << "Computing PLAFs..." << endl;
    ArrayXd plafs = calc_plafs_from_vcf(argv[1]);
    cout << "  PLAFs: " << plafs(0) << "..." << plafs(n_loci - 1) << endl;
}