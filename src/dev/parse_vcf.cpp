#include "htslib/vcf.h"
#include <iostream>
using namespace std;

int main(int argc, char* argv[])
{
    // Exit if incorrect usage
    if (argc != 2) {
        cerr << "Usage: ./a.out <input_file.vcf>" << endl;
        exit(1);
    }

    // Prepare pointers to file, header, and record
    htsFile* vcf_ptr;
    bcf_hdr_t* vcf_hdr_ptr;
    bcf1_t* vcf_rec_ptr = bcf_init(); // for some reason need to initialise

    // Intialise pointers to the file and header
    vcf_ptr = vcf_open(argv[1], "r");
    vcf_hdr_ptr = vcf_hdr_read(vcf_ptr);

    int nsmpl = bcf_hdr_nsamples(vcf_hdr_ptr);
    cout << "Found " << nsmpl << " samples in the header." << endl;

    // Let's iterate over the genotype fields
    int z = 0;
    int ndst = 0;
    int* dst = NULL;
    while (vcf_read(vcf_ptr, vcf_hdr_ptr, vcf_rec_ptr) == 0) {
        // First, we need to unpack
        // if (bcf_unpack(vcf_rec_ptr, BCF_UN_FLT) == 0) {
        //     cout << "Unpacked to 'BCF_UN_FLT'" << endl;
            
            // Get genotypes
            // int load_gt = bcf_get_genotypes(vcf_hdr_ptr, vcf_rec_ptr, &dst, &ndst);
            int load_gt = bcf_get_format_int32(vcf_hdr_ptr, vcf_rec_ptr, "AD", &dst, &ndst);
            if (load_gt > 0) {
                // cout << ndst << "\t";
                // for (int j = 0; j < nsmpl; ++j) {
                //     cout << dst[j] << "\t";
                // }
                cout << "Loaded GT..." << endl;
                cout << "*dst: " << *dst << endl;
                cout << "dst[0]: " << dst[0] << endl;
                cout << "dst[1]: " << dst[1] << endl;
                cout << "dst[2]: " << dst[2] << endl;
                cout << "ndst: " << ndst << endl;
            } else {
                cout << "Failed to load GT, instead got: " << load_gt << endl;
                // We don't simulate genotypes!
            }
        // }

        // Stop after three samples
        ++z;
        if (z > 3) {
            break;
        }

    }
    // free(dst[0]);
    // free(dst);

    // Iterate over records
    // int z = 0;
    // while(vcf_read(vcf_ptr, vcf_hdr_ptr, vcf_rec_ptr) == 0) {
    //     cout << "======================================================================" << endl;
    //     cout << "VARIANT NUMBER: " << z << endl;
    //     cout << "----------------------------------------------------------------------" << endl;
    //     cout << "Direct components of `bcf1_t`:" << endl;
    //     cout << "  ref: " << bcf_hdr_id2name(vcf_hdr_ptr, vcf_rec_ptr->rid) << endl;
    //     cout << "  ->rid: " << vcf_rec_ptr->rid << endl;
    //     cout << "  ->rlen: " << vcf_rec_ptr->rlen << endl;
    //     cout << "  ->pos: " << vcf_rec_ptr->pos << endl;
    //     cout << "  ->n_allele: " << vcf_rec_ptr->n_allele << endl;
    //     cout << "  ->n_sample: " << vcf_rec_ptr->n_sample << endl;
    //     cout << "  ->n_info: " << vcf_rec_ptr->n_info << endl;
    //     cout << "  ->n_fmt: " << vcf_rec_ptr->n_fmt << endl;

    //     cout << "Unpacking until level `BCF_UN_FLT`..." << endl;
    //     // Here, you have to set how far you want to unpack, e.g.
    //     // BCF_UN_FLT  2       // up to FILTER
    //     if (bcf_unpack(vcf_rec_ptr, BCF_UN_FLT) == 0) {
    //         cout << "  Unpacked!" << endl;
    //         cout << "  ->d.als: " << vcf_rec_ptr->d.als << endl;
    //         cout << "  ->d.n_var: " << vcf_rec_ptr->d.n_var << endl;
    //         cout << "  ->d.id: " << vcf_rec_ptr->d.id << endl;
    //         cout << "  ->d.alleles[0]: " << vcf_rec_ptr->d.allele[0] << endl;
    //         cout << "  ->d.alleles[1]: " << vcf_rec_ptr->d.allele[1] << endl;
    //         cout << "  ->d.fmt: " << vcf_rec_ptr->d.fmt << endl;
    //     }

    //     ++z;
    //     if (z > 10) {
    //         break;
    //     }
    // }

    // Relase the memory
    bcf_hdr_destroy(vcf_hdr_ptr);
    bcf_destroy(vcf_rec_ptr);
    vcf_close(vcf_ptr);
}
