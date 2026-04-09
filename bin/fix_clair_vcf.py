#!/usr/bin/env python3
"""
Script for adapting ClairS-TO VCF for viral minor variants.
Restores the quality scores from the sample field to the quality field for failed variants (automatically set to zero).
Removes irrelevant filters (VariantCluster, MultiHap, NoAncestry) from the FILTER column.
Applies LowQual filter if GQ is below 5.
"""

import argparse
import subprocess

# Parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Replace QUAL with GQ from sample field and remove irrelevant filters from VCF."
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Path to the input VCF file."
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Path to the output VCF file."
    )
    return parser.parse_args()

def main():
    args = parse_arguments()

    with open(args.input, "r") as vcf_in, open(args.output, "w") as vcf_out:
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
                filter for filter in filter_col.split(";") 
                if filter not in irrelevant_filters
            ]

            # Add LowQual back in if GQ is below 5 (it gets added if any filter is triggered regardless if the GQ)
            if float(gq_value) < 5:
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
    subprocess.run(["bcftools", "view", args.output, "--output-type", "z", "--output-file", compressed_vcf], check=True)

    # Index it
    subprocess.run(
    ["tabix", "-p", "vcf", compressed_vcf],
    check=True
    )

    print(f"Compressed VCF saved as {compressed_vcf}")

if __name__ == "__main__":
    main()