# Tests for processing of VCF files using vcfR
# JHendry, 2022/11/26

test_that("Converting AD from chr to int", {
    # Separation doesn't matter
    expect_equal(convert_ad_to_int("0, 0"), c(0, 0))
    expect_equal(convert_ad_to_int("0,0"), c(0, 0))
    expect_equal(convert_ad_to_int("  0,   0   "), c(0, 0))
    
    # Various possible hard-coded values
    expect_equal(convert_ad_to_int("0, 53"), c(0, 53))
    expect_equal(convert_ad_to_int("104, 1"), c(104, 1))
    expect_equal(convert_ad_to_int("2348, 0"), c(2348, 0))
    expect_equal(convert_ad_to_int("10, 10"), c(10, 10))
    expect_equal(convert_ad_to_int("99999, 99999"), c(99999, 99999))

    # Possible errors
    expect_error(convert_ad_to_int("0/53"))
    expect_error(convert_ad_to_int("0, 53, 53"))
    expect_error(convert_ad_to_int("0"))
    expect_error(convert_ad_to_int(c(32, 23)))

    # Higher allelic cases
    expect_equal(convert_ad_to_int("0, 53, 53", 3), c(0, 53, 53))
    expect_equal(convert_ad_to_int("0, 102, 10, 5 ", 4), c(0, 102, 10, 5))
    expect_equal(convert_ad_to_int("0, 102, 10, 5, 6", 5), c(0, 102, 10, 5, 6))
})


test_that("Loading and processing VCF data", {
    # Load an example VCF
    eg_vcf = "../../example_data/simulated_infections.DRCongo.K02.vcf"
    expect_true(file.exists(eg_vcf))
    sample_data <- load_sample_vcf_data(eg_vcf, target_sample=1)

    # Structure of list
    expected_names <- c("chroms", "pos", "refs", "alts", "plafs", "wsafs")
    expect_equal(names(sample_data), expected_names)

    # Sanity check length of list elements
    lengths <- unlist(lapply(sample_data, length))
    expect_true(sum(lengths) > 0)
    expect_equal(min(lengths), max(lengths))

    # TODO:
    # - Include checks of the data types
    # - Include behaviours from *broken* VCF files

})
