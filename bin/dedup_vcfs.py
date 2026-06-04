#!/usr/bin/env python3
"""
Script parses a consensus level VCF and a ClairS-TO VCF to remove redundant variants.
Input: a consensus VCF and a ClairS-TO minor variant VCF from the same sample.
Output: VCF containing variants unique to ClairS-TO (i.e. removes redundant variants that are present in the consensus VCF).
Splits multi-nucleotide variants into single-nucleotide variants in the consensus VCF (necessary due to Medaka).
Uses bcftools isec function to output unique variants (must be installed).
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
        "Script parses a Consensus VCF and a ClairS-TO VCF to remove redundant variants."
        "Input: a consensus VCF and a ClairS-TO minor variant VCF."
        "Output: VCF containing minor variants unique to ClairS-TO."
    ))
    parser.add_argument("--consensus-vcf", required=True, help="Path to the Consensus VCF file.")
    parser.add_argument("--clairSTO-vcf", required=True, help="Path to the ClairS-TO VCF file.")
    return parser

def split_mnv_to_snvs(record: pysam.VariantRecord) -> list:
    """
    Splits a multi-nucleotide variant (MNV) VCF entry into single-nucleotide variants (SNVs).
    Expects input as records of a VCF supplied from pysam.VariantFile.
    Returns VCF records in list.
    """
    snvs = []
    ref = record.ref
    alt = record.alts[0]

    # check for multiallelic sites (using normalized vcf should prevent this!)
    if len(record.alts) > 1:
        print(f"Skipping multiallelic variant at {record.chrom}:{record.pos}. Run BCFtools norm to decompose these!")
        return [record]

    # Redundant check for indels
    if len(ref) != len(alt):
        print(f"Skipping non-substitution variant at {record.chrom}:{record.pos}.")
        return [record]

    # Iterate over bases of the variant and print as individual records
    for i, (ref_base, alt_base) in enumerate(zip(ref, alt)):
        snv = record.copy()
        snv.pos = record.pos + i
        snv.ref = ref_base
        snv.alts = (alt_base,)
        snvs.append(snv)

    return snvs

def process_vcf(consensus_vcf: str, clairSTO_vcf: str) -> None:
    """
    Processes the consensus VCF to split MNVs into SNVs.
    Deduplicates VCFs using bcftools isec
    """

    base_name = os.path.basename(consensus_vcf)
    if base_name.endswith('.vcf.gz'):
        base_name = base_name[:-7]
    elif base_name.endswith('.vcf'):
        base_name = base_name[:-4]

    # Intermediate files
    split_vcf = f"{base_name}.split.vcf"
    compressed_vcf = f"{base_name}.split.vcf.gz"
    sorted_vcf = f"{base_name}.sorted.vcf.gz"

    # Split MNVs into SNVs in the consensus vcf
    with pysam.VariantFile(consensus_vcf, 'r') as vcf_reader, pysam.VariantFile(split_vcf, 'w', header=vcf_reader.header) as vcf_writer:
        for record in vcf_reader:
            # find and split MNVs but make sure they aren't indels or multiallelic, print everything else as is
            if len(record.alts) == 1 and len(record.ref) == len(record.alts[0]) and len(record.ref) > 1:
                snvs = split_mnv_to_snvs(record)
                for snv in snvs:
                    vcf_writer.write(snv)
            else:
                vcf_writer.write(record)
    print(f"Split MNVs written to {split_vcf}")

    # Compress, sort and index the consensus VCF (needed for isec)
    subprocess.run(["bcftools", "view", split_vcf, "--output-type", "z", "--output-file", compressed_vcf], check=True)
    subprocess.run(["bcftools", "sort", compressed_vcf, "--output-type", "z", "--output-file", sorted_vcf], check=True)
    subprocess.run(["tabix", "-p", "vcf", sorted_vcf], check=True)
    print(f"Compressed and sorted VCF written to {sorted_vcf}")

    # Deduplicate using bcftools and the clairSTO vcf
    dedup_dir = "dedup_vcfs"
    subprocess.run(["bcftools", "isec", "-p", dedup_dir, clairSTO_vcf, sorted_vcf], check=True)
    unique_vcf = os.path.join(dedup_dir, "0000.vcf") # this vcf has variants unique to ClairS-TO
    print(f"Deduplicated ClairS-TO VCF written to {unique_vcf}")

if __name__ == "__main__":
    parser = init_parser()
    args = parser.parse_args()
    process_vcf(args.consensus_vcf, args.clairSTO_vcf)
