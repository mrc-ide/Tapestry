

// Other installs
#include <iostream>

// Can we manage to include HTSLIB?
// #include "/usr/local/include/htslib/vcf.h"
#include "htslib/vcf.h"
using namespace std;

int main(int argc, char *argv[])
{
    htsFile *fp;
    bcf_hdr_t *hdr;
    bcf_idpair_t *ctg;
    int i;
    if (argc == 1) {
        fprintf(stderr, "Usage: print-ctg <in.vcf>\n");
        return 1;
    }
    fp = vcf_open(argv[1], "r");
    hdr = vcf_hdr_read(fp);
    ctg = hdr->id[BCF_DT_CTG];
    for (i = 0; i < hdr->n[BCF_DT_CTG]; ++i) {
        cout << ctg[i].key << endl;
        cout << ctg[i].val->info[0] << endl;
    }
    //printf("%s\t%f\n", ctg[i].key, ctg[i].val->info[0]);
    bcf_hdr_destroy(hdr);
    vcf_close(fp);
    return 0;
}