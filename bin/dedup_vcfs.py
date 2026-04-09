#!/usr/bin/env python3
"""
Script parses a consensus level VCF and a ClairS-TO VCF to remove redundant variants.
Input: a consensus VCF and a ClairS-TO minor variant VCF from the same sample.
Output: VCF containing variants unique to ClairS-TO (i.e. removes redundant variants that are present in the consensus VCF).
Splits multi-nucleotide variants into single-nucleotide variants in the consensus VCF (necessary due to Medaka).
Uses bcftool's isec function to output unique variants (must be installed).
"""

import argparse
import pysam
import subprocess
import os

def init_parser() -> argparse.ArgumentParser:
    """
    Command-line argument parsing.
    """
    parser = argparse.ArgumentParser(description=(
        "Scrip parses a Medaka VCF and a ClairS-TO VCF to remove redundant variants."
        "Input: a Medaka consensus VCF and a ClairS-TO minor variant VCF."
        "Output: VCF containing minor variants unique to ClairS-TO."
    ))
    parser.add_argument("--medaka-vcf", required=True, help="Path to the Medaka VCF file.")
    parser.add_argument("--clairs-vcf", required=True, help="Path to the Clairs-to VCF file.")
    return parser

def split_mnv_to_snvs(record):
    """
    Splits a multi-nucleotide variant (MNV) VCF entry into single-nucleotide variants (SNVs).
    Expects input as records of a Medaka VCF supplied from pysam.VariantFile.
    Returns VCFs records one by one.
    """
    snvs = []
    ref = record.ref
    alt = record.alts[0]  
    
    # Redundant check for indels
    if len(ref) != len(alt):
        print(f"Skipping non-substitution variant at {record.chrom}:{record.pos} (REF={ref}, ALT={alt}).")
        return snvs
    
    # Iterate over bases of the variant and print as individual records 
        # (zip produces pairs of ref-alt and enumerate indexes and returns tuples)
    for i, (ref_base, alt_base) in enumerate(zip(ref, alt)):
        snv = record.copy()
        snv.pos = record.pos + i
        snv.ref = ref_base
        snv.alts = (alt_base,)
        snvs.append(snv)

    return snvs

def process_vcf(medaka_vcf, clairs_vcf):
    """
    Processes the Medaka VCF to split MNVs into SNVs.
    Deduplicates VCFs using bcftools isec
    """

    base_name = os.path.basename(medaka_vcf)
    if base_name.endswith('.vcf.gz'):
        base_name = base_name[:-7]
    elif base_name.endswith('.vcf'):
        base_name = base_name[:-4]

    # Intermediate files
    split_vcf = f"{base_name}.split.vcf"
    compressed_vcf = f"{base_name}.split.vcf.gz"
    sorted_vcf = f"{base_name}.sorted.vcf.gz"

    # Split MNVs into SNVs in the medaka vcf
    with pysam.VariantFile(medaka_vcf, 'r') as vcf_reader, pysam.VariantFile(split_vcf, 'w', header=vcf_reader.header) as vcf_writer:
        for record in vcf_reader:
            # find and split MNVs but make sure they aren't indels, print everything else
            if len(record.ref) == len(record.alts[0]) and len(record.ref) > 1:
                snvs = split_mnv_to_snvs(record)
                for snv in snvs:
                    vcf_writer.write(snv)
            else:
                vcf_writer.write(record)
    print(f"Split MNVs written to {split_vcf}")

    # Compress, sort and index the consensus VCF (needed for ised)
    subprocess.run(["bcftools", "view", split_vcf, "--output-type", "z", "--output-file", compressed_vcf], check=True)
    subprocess.run(["bcftools", "sort", compressed_vcf, "--output-type", "z", "--output-file", sorted_vcf], check=True)
    subprocess.run(["tabix", "-p", "vcf", sorted_vcf], check=True)
    print(f"Compressed and sorted VCF written to {sorted_vcf}")

    # Deduplicate using bcftools and the clairs-to vcf
    dedup_dir = "dedup_vcfs"
    subprocess.run(["bcftools", "isec", "-p", dedup_dir, clairs_vcf, sorted_vcf], check=True)
    unique_vcf = os.path.join(dedup_dir, "0000.vcf") # this vcf has variants unique to ClairS-TO
    print(f"Deduplicated Clairs-to VCF written to {unique_vcf}")

    # Cleanup temp files
    os.remove(split_vcf)
    os.remove(compressed_vcf)

if __name__ == "__main__":
    parser = init_parser()
    args = parser.parse_args()
    process_vcf(args.medaka_vcf, args.clairs_vcf)