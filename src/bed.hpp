#pragma once
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>


struct BEDRecord
{
private:
    void check_valid();

public:
    // Members
    std::string chrom;
    int start;
    int end;
    std::string name;

    // Constructors
    BEDRecord();
    BEDRecord(std::string chrom, int start, int end);
    BEDRecord(std::string chrom, int start, int end, std::string name);

    friend std::ostream& operator<<(std::ostream& os, const BEDRecord& bed)
    {
        os << bed.chrom << "\t" << bed.start << "\t" << bed.end;
        if (bed.name.length() > 0) {
            os << "\t" << bed.name;
        }
        os << "\n";

        return os;
    }
};


void write_bed_records(
    const std::string& output_bed,
    const std::vector<BEDRecord>& bed_records
);

