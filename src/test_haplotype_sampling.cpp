#include <vector>
#include <numeric>
#include <cmath>
#include <testthat.h>
#include "combinatorics.h"
#include "haplotype_sampling.h"
using namespace std;



/**
 * To test:
 * - Size of array X
 * - Probabilities always sum to 1.0 X (working, but assertion fails)
 * - Probabilities are always in the expected p^1:N
*/


context("Test haplotype sampling probabilities")
{
    // Create range of COI values over which to test
    int MAX_K = 10;
    vector<int> Ks(MAX_K - 1);
    std::iota(Ks.begin(), Ks.end(), 1);

    test_that("Test size of powerset matrix")
    {

        double p = 0.5;

        for (int k : Ks)
        {
            // Create strain indices
            vector<int> strains(k);
            std::iota(strains.begin(), strains.end(), 0);

            // Haplotype and IBD configurations
            vector<vector<int>> hap_configs;
            hap_configs = create_powerset(k);

            vector<vector<vector<int>>> ibd_configs;
            ibd_configs = create_all_partitions(strains);

            // Create haplotype sampling matrix
            vector<vector<double>> sampling_probs;
            sampling_probs = calc_sampling_probs(p, ibd_configs, hap_configs);

            // Check size
            // DO I WANT THIS ORIENNTATION?
            // HAP CONFIGS ARE ROWS
            expect_true(sampling_probs.size() == hap_configs.size());
            expect_true(sampling_probs[0].size() == ibd_configs.size());
        }
    }

    test_that("Given an IBD configuration, haplotype sampling probabilities should sum to 1.")
    {

        // A bunch of allele frequencies to try
        //vector<double> ps = {0, 0.01, 0.1, 0.2, 0.5, 0.85, 0.97, 0.4444, 1.0};
        vector<double> ps = {0.1, 0.25, 0.5, 0.75, 0.9};
        for (double p : ps)
        {
            for (int k : Ks)
            {
                // Create strain indices
                vector<int> strains(k);
                std::iota(strains.begin(), strains.end(), 0);

                // Haplotype and IBD configurations
                vector<vector<int>> hap_configs;
                hap_configs = create_powerset(k);

                vector<vector<vector<int>>> ibd_configs;
                ibd_configs = create_all_partitions(strains);

                // Create haplotype sampling matrix
                vector<vector<double>> sampling_probs;
                sampling_probs = calc_sampling_probs(p, ibd_configs, hap_configs);

                for (int j = 0; j < ibd_configs.size(); ++j)
                {
                    double total_prob = 0.0;
                    for (int i = 0; i < hap_configs.size(); ++i)
                    {
                        total_prob += sampling_probs[i][j];
                    }
                    expect_true(total_prob == 1.0);
                }
            }
        }
    }
}