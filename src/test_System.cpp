#include <testthat.h>
#include "System.h"



/**
 * Goal is to test that SystemVCF is instantiated as expected
 * In effect, that the data type transformations work well.
 * But for this, I need to start in R to load the VCF and then
 * create the SystemVCF object; how to do this?
*/
context("Test loading of VCF data into System") {
    test_that("Addition") {
        int a = 1;
        int b = 2;
        expect_equal(a + b == 3);
    }
}