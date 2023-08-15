#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include "data.hpp"
#include "ibd.hpp"
#include "libs/eigen-3.4.0/Dense"
using Eigen::VectorXd;
using Eigen::MatrixBase; // Base class of all matrix types; for templating


void write_proportions(
    const string& output_csv, 
    const VectorXd ws, 
    double epsilon = 0.0001);


template<typename Derived>
void write_data_with_annotation(
    std::string output_csv,
    const VCFData& data,
    const MatrixBase<Derived>& annotation_matrix,
    const vector<string>& annotation_columns
)
{
    // Check shapes match
    if (annotation_matrix.rows() != data.n_sites) {
        throw std::invalid_argument("Annotation matrix has incorrect number of rows.");
    }
    
    // Open the output file
    std::ofstream csv_file(output_csv);
    if (!csv_file.is_open()) {
        throw std::invalid_argument("Could not open output file.");
    }

    // Write column names
    csv_file << "chrom,pos,refs,alts,wsafs,plafs";
    for (const string& col_name : annotation_columns) {
        csv_file << "," << col_name;
    }
    csv_file << "\n";

    // Format outputs
    const static Eigen::IOFormat CSVFormat(
        5,
        Eigen::DontAlignCols,
        ", ", 
        "\n",
        "",
        "\n"
        );

    // Write data
    for (int i = 0; i < data.n_sites; ++i) {
        csv_file << data.chrom_names[i] << ",";
        csv_file << data.pos(i) << ",";
        csv_file << data.refs(i) << ",";
        csv_file << data.alts(i) << ",";
        csv_file << data.wsafs(i) << ",";
        csv_file << data.plafs(i) << ",";
        csv_file << annotation_matrix.row(i).format(CSVFormat);
    }

    csv_file.close();
}

