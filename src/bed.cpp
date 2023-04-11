#include "bed.hpp"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>


BEDRecord::BEDRecord() {}

BEDRecord::BEDRecord(std::string chrom, int start, int end)
    : chrom(chrom), 
    start(start), 
    end(end)
{
    check_valid();
}

BEDRecord::BEDRecord(std::string chrom, int start, int end, std::string name)
    : chrom(chrom),
    start(start),
    end(end),
    name(name)
{
    check_valid();
}


void BEDRecord::check_valid()
{
    if (end < start) {
        throw  std::invalid_argument("Start position must be less than end position.");
    }
}


void write_bed_records(
    const std::string& output_bed,
    const std::vector<BEDRecord>& bed_records
)
{
     // Open the output file
    std::ofstream bed_file(output_bed);
    if (!bed_file.is_open()) {
        throw std::invalid_argument("Could not open output file.");
    }

    for (const BEDRecord& bed : bed_records) {
        bed_file << bed;
    }

    bed_file.close();
}