#include <iostream>
#include <string>
#include <vector>
#include "data.hpp"
#include "misc_v15.h"
using namespace std;


Data::Data(
    int n_loci,
    Rcpp::StringVector chroms,
    Rcpp::NumericVector pos,
    Rcpp::NumericVector refs,
    Rcpp::NumericVector alts,
    Rcpp::NumericVector plafs,
    Rcpp::NumericVector wsafs
)
    : n_loci(n_loci),
    chroms(rcpp_to_vector_string(chroms)), 
    pos(rcpp_to_vector_int(pos)), 
    refs(rcpp_to_vector_int(refs)), 
    alts(rcpp_to_vector_int(alts)), 
    plafs(rcpp_to_vector_double(plafs)), 
    wsafs(rcpp_to_vector_double(wsafs))
{
}


void Data::print()
{
    cout << "Data" << endl;
    cout << "  Chromosomes: " << chroms[0] << " ... " << chroms[n_loci - 1] << endl;
    cout << "  Positions: " << pos[0] << " ... " << pos[n_loci - 1] << endl;
    cout << "  REF counts: " << refs[0] << " ... " << refs[n_loci - 1] << endl;
    cout << "  ALT counts: " << alts[0] << " ... " << alts[n_loci - 1] << endl;
    cout << "  PLAFs: " << plafs[0] << " ... " << plafs[n_loci - 1] << endl;
    cout << "  WSAFs " << wsafs[0] << " ... " << wsafs[n_loci - 1] << endl;
}

  

   