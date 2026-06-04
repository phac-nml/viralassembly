#!/usr/bin/env python3
"""
Script for adapting ClairS-TO VCF for viral minor variants.
Restores the quality scores from the sample field to the quality field for failed variants (automatically set to zero when a filter is failed).
Removes irrelevant filters (VariantCluster, MultiHap, NoAncestry) from the FILTER column.
Applies LowQual filter if GQ is below 5 or <quality>.
Reliant on bcftools for compression and indexing.

"""

import argparse
import pysam
import gzip
import os

# Parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Replace QUAL with GQ from sample field and remove filters not applicable to viruses from VCF."
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Path to the input VCF file."
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Path to the output VCF file."
    )
    parser.add_argument(
        "-q", "--quality", required=False, default=5, type=int, help="Minimum quality score to pass a variant."
    )
    return parser.parse_args()

def fix_non_ascii(input_vcf: str) -> str:
    """
    Reads in vcf and replaces non-ASCII with ?
    Returns path to a temp clean vcf file
    """
    # Temp vcf
    if input_vcf.endswith('.vcf.gz'):
        vcf_fix = input_vcf.replace(".vcf.gz", "_clean.vcf")
    else:
        vcf_fix = input_vcf.replace(".vcf", "_clean.vcf")
    
    # In case input is zipped
    zip = input_vcf.endswith('.gz')
    open_func = gzip.open if zip else open
    read_mode = 'rt' if zip else 'r'

    # Clean any non-ASCII characters from vcf
    with open_func(input_vcf, read_mode, encoding='utf-8', errors='replace') as vcf_in, \
    open(vcf_fix, 'w', encoding='ascii', errors='replace') as vcf_out:
        for line in vcf_in:
            vcf_out.write(line)

    return vcf_fix

def main() -> None:
    args = parse_arguments()

    irrelevant_filters = {"VariantCluster", "MultiHap", "NoAncestry", "NonSomatic", "Realignment"}

    # Fix non-ASCII characters because ClairS headers
    cleaned_vcf = fix_non_ascii(args.input)

    with pysam.VariantFile(cleaned_vcf, 'r') as vcf_in:

        new_header = vcf_in.header.copy()

        # Remove irrelevent filters (i.e. for euks/diploid)
        for filter_id in irrelevant_filters:
            if filter_id in new_header.filters:
                new_header.filters.remove_header(filter_id)

        output_file = args.output if args.output.endswith('.gz') else f"{args.output}.gz"
        with pysam.VariantFile(output_file, 'wz', header=new_header) as vcf_out:
            for record in vcf_in:
                
                # Grab Qscore from info field
                gq_value = record.samples[0].get("GQ", 0)

                # Return Qscore to Qual field as filter automatically sets to 0
                record.qual = float(gq_value)

                # Remove irrelevant filters
                irrelevant_filters.add("LowQual")  # Add LowQual here because it gets added for any filter pass, add back based on our qual threshold later
                current_filters = set(record.filter)
                remaining_filters = current_filters - irrelevant_filters

                # Add LowQual back if below threshold
                if gq_value < args.quality:
                    remaining_filters.add("LowQual")

                # Reset filters (empty = PASS)
                record.filter.clear()
                if remaining_filters:
                    for f in remaining_filters:
                        record.filter.add(f)
                else:
                    record.filter.add("PASS")
                
                vcf_out.write(record)
    
    # Index the VCF
    pysam.tabix_index(output_file, preset="vcf", force=True)
    
    print(f"Updated VCF saved as {output_file}")


if __name__ == "__main__":
    main()
