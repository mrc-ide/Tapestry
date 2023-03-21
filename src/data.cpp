#include <iostream>
#include <string>
#include <vector>
#include "vcf.hpp"
#include "data.hpp"
using namespace std;


VCFData::VCFData(const string& vcf_path, const string& sample_name)
    : vcf_path(vcf_path),
    sample_name(sample_name),
    n_sites(count_sites_in_vcf(vcf_path)),
    chrom_names(n_sites),
    chroms(ArrayXi::Constant(n_sites, -1)),
    pos(ArrayXi::Constant(n_sites, -1)),
    refs(ArrayXi::Constant(n_sites, -1)),
    alts(ArrayXi::Constant(n_sites, -1)),
    wsafs(ArrayXd::Constant(n_sites, -1.0)),
    plafs(ArrayXd::Constant(n_sites, -1.0))
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
            if (sz != 2) {
                throw std::invalid_argument("Error: all sites must be bi-allelic.");
            }
            refs(site_ix) = vals[0];
            alts(site_ix) = vals[1];
            wsafs(site_ix) = vals[1] / (vals[0] + vals[1] + epsilon);
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

    // Calculate PLAFs
    plafs = calc_plafs_from_vcf(vcf_path);

    // TODO:
    // - Check all the output arrays are sensible (i.e. in correct ranges)
}


void VCFData::print() const
{
    cout << "Loaded VCF data" << endl;
    cout << "  VCF path: " << vcf_path << endl;
    cout << "    No. Sites: " << n_sites << endl;
    cout << "    No. Samples: " << n_samples << endl;
    cout << "  Target sample: " << sample_name << endl;

    // TODO: should make a macro or inline function for these...
    cout << "    Chromosomes: " << chrom_names[0] << " ... " << chrom_names[n_sites - 1] << endl;
    cout << "    Positions: " << pos(0) << " ... " << pos(n_sites - 1) << endl;
    cout << "    REFs: " << refs(0) << " ... " << refs(n_sites - 1) << endl;
    cout << "    ALTs: " << alts(0) << " ... " << alts(n_sites - 1) << endl;
    cout << "    WSAFs: " << wsafs(0) << " ... " << wsafs(n_sites - 1) << endl;
    cout << "    PLAFs: " << plafs(0) << " ... " << plafs(n_sites - 1) << endl;
}

