
#include "System.h"
#include "misc_v15.h"
#include "typedefs.hpp"
#include "combinatorics.h"
#include "haplotype_sampling.h"

using namespace std;




//--------------------------------------------------------------------------------
// Public member functions
//
//--------------------------------------------------------------------------------



SystemVCF::SystemVCF(
    int K, 
    int n_loci,
    Rcpp::StringVector chroms,
    Rcpp::NumericVector pos,
    Rcpp::NumericVector refs,
    Rcpp::NumericVector alts,
    Rcpp::NumericVector plafs,
    Rcpp::NumericVector wsafs
  )
    : K(K), 
    e_0(0.01),
    e_1(0.05),
    v(500),
    n_loci(n_loci), 
    chroms(rcpp_to_vector_string(chroms)), 
    pos(rcpp_to_vector_int(pos)), 
    refs(rcpp_to_vector_int(refs)), 
    alts(rcpp_to_vector_int(alts)), 
    plafs(rcpp_to_vector_double(plafs)), 
    wsafs(rcpp_to_vector_double(wsafs))
{
    print();
}


SystemVCF::SystemVCF(
    int K, 
    int n_loci,
    double e_0,
    double e_1,
    double v,
    Rcpp::StringVector chroms,
    Rcpp::NumericVector pos,
    Rcpp::NumericVector refs,
    Rcpp::NumericVector alts,
    Rcpp::NumericVector plafs,
    Rcpp::NumericVector wsafs
  )
    : K(K), 
    e_0(e_0),
    e_1(e_1),
    v(v),
    n_loci(n_loci), 
    chroms(rcpp_to_vector_string(chroms)), 
    pos(rcpp_to_vector_int(pos)), 
    refs(rcpp_to_vector_int(refs)), 
    alts(rcpp_to_vector_int(alts)), 
    plafs(rcpp_to_vector_double(plafs)), 
    wsafs(rcpp_to_vector_double(wsafs))
{
    print();
}


void SystemVCF::print()
{
    cout << "Data:" << endl;
    cout << "  Chromosomes: " << chroms[0] << " ... " << chroms[n_loci - 1] << endl;
    cout << "  Positions: " << pos[0] << " ... " << pos[n_loci - 1] << endl;
    cout << "  REF counts: " << refs[0] << " ... " << refs[n_loci - 1] << endl;
    cout << "  ALT counts: " << alts[0] << " ... " << alts[n_loci - 1] << endl;
    cout << "  PLAFs: " << plafs[0] << " ... " << plafs[n_loci - 1] << endl;
    cout << "  WSAFs " << wsafs[0] << " ... " << wsafs[n_loci - 1] << endl;
    cout << "Parameters:" << endl;
    cout << "  K: " << K << endl;
    cout << "  e_0: " << e_0 << endl;
    cout << "  e_1: " << e_1 << endl;
    cout << "  v: " << v << endl;
    cout << "Done." << endl;
}




//--------------------------------------------------------------------------------
// Private member functions
//
//--------------------------------------------------------------------------------



