#include <cmath>
#include <numeric>
#include "System.h"
#include "misc_v15.h"
#include "typedefs.hpp"
#include "combinatorics.h"
#include "haplotype_sampling.h"
#include "betabin.hpp"
using namespace std;


//--------------------------------------------------------------------------------
// System : Public member functions
//
//--------------------------------------------------------------------------------


System::System(Data& data, Parameters& params) :
  data(data), 
  params(params),
  strains(params.K),
  hap_configs(boost::extents[pow(2, params.K)][params.K])
{
  print();
}

void System::print()
{
  data.print();
  params.print();
}

void System::precompute_arrays()
{
  create_strains();
  create_hap_configs();
  create_ibd_configs();
  create_hap_sampling_probs();
  create_betabin_array();
}


//--------------------------------------------------------------------------------
// Private member functions
//
//--------------------------------------------------------------------------------


void System::create_strains()
{
  iota(strains.begin(), strains.end(), 0);
}


void System::create_hap_configs()
{
  hap_configs = create_powerset(params.K);
}


void System::create_ibd_configs()
{
  ibd_configs = create_all_partitions(strains);
}


void System::create_hap_sampling_probs()
{
  // We produce one haplotype sampling array per PLAF
  int n_plafs = data.plafs.size();
  this->hap_sampling_probs.resize(boost::extents[n_plafs][hap_configs.size()][ibd_configs.size()]);
  
  // Compute and assign here
  for (int i = 0; i < n_plafs; ++i) {
      this->hap_sampling_probs[i] = calc_sampling_probs(data.plafs[i], ibd_configs, hap_configs);
  }
}

// TODO:
// - This is bad in a lot of ways
// - First, I should inject
// - Second, I don't want to instantiate and then reassign
// - Really, this should come in at the constructor
void System::create_betabin_array()
{

  int n_pi_bins = 100;
  BetabinomialArrayByHash bbarray(
    data.n_loci, 
    data.refs, 
    data.alts, 
    n_pi_bins, 
    params.e_0, 
    params.e_1, 
    params.v
  );

  cout << "Computing betabinomial lookup...";
  bbarray.compute_array();
  // this->betabin_array = bbarray;
  cout << " Done." << endl;
}
