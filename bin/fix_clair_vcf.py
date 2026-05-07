#!/usr/bin/env python3
"""
Script for adapting ClairS-TO VCF for viral minor variants.
Restores the quality scores from the sample field to the quality field for failed variants (automatically set to zero when a filter is failed).
Removes irrelevant filters (VariantCluster, MultiHap, NoAncestry) from the FILTER column.
Applies LowQual filter if GQ is below 5 or <quality>.
Reliant on bcftools for compression and indexing.

"""

import argparse
import subprocess

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

def main() -> None:
    args = parse_arguments()

    with open(args.input, "r", encoding='utf-8', errors='replace') as vcf_in, \
         open(args.output, "w", encoding='ascii', errors='replace') as vcf_out: # Ensure ASCII encoding for output otherwise downstream issues parsing the vcf
        for line in vcf_in:
            # Write header lines, but skip irrelevant filter definitions
            if line.startswith("#"):
                if any(line.startswith(f"##FILTER=<ID={warn}") for warn in ["VariantCluster", "MultiHap", "NoAncestry", "NonSomatic", "Realignment"]):
                    continue
                vcf_out.write(line)
                continue

            # Process variants, strip any white space and tab is seperater
            fields = line.strip().split("\t")

            # Replace QUAL with GQ
            sample_column = fields[9]
            # Info field split by : and gq is second (GT:GQ:DP:AF:AD:AU:CU:GU:TU)
            gq_value = sample_column.split(":")[1]
            fields[5] = gq_value

            # Remove irrelevant filters from FILTER column
            filter_col = fields[6]
            irrelevant_filters = {"VariantCluster", "MultiHap", "NoAncestry", "LowQual"}
            filters = [
                filt for filt in filter_col.split(";")
                if filt not in irrelevant_filters
            ]

            # Add LowQual back in if GQ is below args.quality (Default 5) (it gets added if any filter is triggered regardless of the GQ)
            if float(gq_value) < args.quality:
                if "LowQual" not in filters:
                    filters.append("LowQual")

            # If no filters left set to PASS
            if not filters:
                fields[6] = "PASS"
            else:
                fields[6] = ";".join(filters)

            # Write the modified fields to the output
            vcf_out.write("\t".join(fields) + "\n")

    print(f"Updated VCF saved as {args.output}")

    # Compress the output VCF with bcftools
    compressed_vcf = f"{args.output}.gz"
    result = subprocess.run(
        ["bcftools", "view", args.output, "--output-type", "z", "--output-file", compressed_vcf],
        check=True,
        capture_output=True,
        text=True
    )
    if result.stderr:
        print(f"bcftools warning: {result.stderr}")

    # Index it
    result = subprocess.run(
        ["tabix", "-f", "-p", "vcf", compressed_vcf],
        check=True,
        capture_output=True,
        text=True
    )
    if result.stderr:
        print(f"tabix warning: {result.stderr}")

    print(f"Compressed VCF saved as {compressed_vcf}")

if __name__ == "__main__":
    main()
    