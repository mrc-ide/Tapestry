#pragma once

// Bell numbers one to ten; from https://oeis.org/A000110
const int BELL_NUMBERS[11] = {
        1,
        1,
        2,
        5,
        15,
        52,
        203,
        877,
        4140,
        21147,
        115975
};

// Value assigned if allelic depth (AD) data is missing
const int MISSING_AD_VALUE = -1;

const int DEFAULT_INT = -9999;
const float DEFAULT_FLOAT = -9999.0;
