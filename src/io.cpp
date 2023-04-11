#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include "data.hpp"
#include "ibd.hpp"
#include "io.hpp"
#include "libs/eigen-3.4.0/Dense"
using Eigen::VectorXd;

/* Write a single set of proportions to a CSV
*  Pre: `prop_csv` should exist; `ws` should sum to 1.
*  Post: A CSV is written to `prop_csv`.
*/
void write_proportions(const string& output_csv, const VectorXd ws, double epsilon)
{
    
    // Check we are writing something that looks like proportions
    if (fabs(ws.sum() - 1.0) > epsilon) {
        throw std::invalid_argument("Proportions do not sum to 1.");
    }
    
    // Create file, check open
    ofstream csv_file(output_csv);
    if (!csv_file.is_open()) {
        throw std::invalid_argument("Could not open file.");
    }

    // Write
    csv_file << "strain,prop\n";
    for (int i = 0; i < ws.size(); ++i) {
        csv_file << i << "," << ws[i] << "\n";
    }

    // Close
    csv_file.close();
}
