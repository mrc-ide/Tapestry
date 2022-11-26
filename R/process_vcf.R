# Start by writing psuedo-code
library(vcfR)


#------------------------------------------------
#' @title Convert string AD representation to numeric
#'
#' @description vcfR's extract.gt function loads genotype values,
#' include allelic depth (AD) information as strings, e.g. '0,53'.
#' Convert these to a numeric vector c(0, 53).
#'
#' @param ad_as_char; string; character representation
#' @return ad_as_char; vector, numeric; numeric representation
convert_ad_to_int <- function(ad_as_char) {
  
  # Handle case where missing
  if (is.na(ad_as_char)) {
    return(NA)
  }
  
  ad_as_ints <- as.numeric(unlist(strsplit(ad_as_char, ",")))

  return(ad_as_ints)
}


#------------------------------------------------
#' @title Load VCF data for a specific sample
#'
#' @description Load VCF data required to run Tapestry for a specific sample.
#'
#'TODO: 
#' - How to best compute PLAF
#' - How to handle sites with zero depth?
#' - How to handle NAs in AD?
#'
#' @param vcf_path; character; path to VCF file.
#' @param target_sample; character or int; sample for which to load data, 
#' specified by either the sample name or index. 
#' @return _; list; VCF data required to run Tapestry, including
#' for every SNP, the chromosome, position, REF read counts, ALT read counts
#' and within-sample allele frequency.
#' @export
load_sample_vcf_data <- function(vcf_path, target_sample) {

    # Load
    vcf <- read.vcfR(vcf_path)

    # Extract chromosome and position
    chroms <- vcf@fix[, "CHROM"]
    pos <- as.numeric(vcf@fix[, "POS"])

    # Load allelic depths for all samples
    ads_chr <- extract.gt(vcf, element="AD")

    # Convert to numeric
    # NB: vcfR's as.numeric=TRUE flag does not work
    ads <- apply(ads_chr, c(1, 2), convert_ad_to_int)

    # Compute PLAFs
    plafs <- rowSums(ads[2, , ])/rowSums(colSums(ads))

    # Extract sample REF and ALT
    refs <- ads[1, , target_sample]
    alts <- ads[2, , target_sample]
    depth <- refs + alts

    # Compute WSAF
    wsafs <- alts / (depth + 1)

    return(list(
        "chroms" = chroms,
        "pos" = pos,
        "refs" = refs,
        "alts" = alts,
        "plafs" = plafs,
        "wsafs" = wsafs
    ))
}








