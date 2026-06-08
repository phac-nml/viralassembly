#!/usr/bin/env python3
'''
Filter VCF Variants originally from https://github.com/artic-network/fieldbioinformatics/blob/master/artic/vcf_filter.py

Notable Changes:
  - Added a custom parameter to adjust the QUAL threshold with it set at 8 by default.
  - Added a filter for RefCall bases.
'''

from cyvcf2 import VCF, Writer
from collections import defaultdict

def in_frame(v):
    if len(v.ALT) > 1:
       print ("This code does not support multiple genotypes!")
       raise SystemExit
    ref = v.REF
    alt = v.ALT[0]
    bases = len(alt) - len(ref)
    if not bases:
       return True
    if bases % 3 == 0:
       return True
    return False


class NanopolishFilter:
    def __init__(self, no_frameshifts):
        self.no_frameshifts = no_frameshifts
        pass

    def check_filter(self, v):
        total_reads = float(v.INFO['TotalReads'])
        qual = v.QUAL
        strandbias = float(v.INFO['StrandFisherTest'])

        if qual / total_reads < 3:
            return False

        if self.no_frameshifts and not in_frame(v):
            return False

        if v.is_indel:
            strand_fraction_by_strand = v.INFO['SupportFractionByStrand']
            if float(strand_fraction_by_strand[0]) < 0.5:
                return False

            if float(strand_fraction_by_strand[1]) < 0.5:
                return False

        if total_reads < 20:
            return False

        return True


class MedakaFilter:
    def __init__(self, no_frameshifts):
        self.no_frameshifts = no_frameshifts

    def check_filter(self, v):
        depth = v.INFO['DP']
        if depth < 20:
            return False

        if self.no_frameshifts and not in_frame(v):
            return False

        if v.num_het:
            return False
        return True


class Clair3Filter:
    def __init__(self, no_frameshifts, min_depth, min_variant_qual, min_frameshift_qual, min_allele_freq):
        self.no_frameshifts = no_frameshifts
        self.min_depth = min_depth
        self.min_variant_qual = min_variant_qual
        self.min_frameshift_qual = min_frameshift_qual
        self.min_allele_freq = min_allele_freq

    def check_filter(self, v):
        qual = v.QUAL

        # Filter LowQual variants
        if qual is None:
            return False

        # Filter out low allele frequency variants
        try:
            allele_freq = v.format("AF")[0][0]
        except Exception:
            print(
                f"ERROR: Could not find AF for variant at {v.CHROM}:{v.POS}, cannot filter on allele frequency"
            )
            raise SystemExit(1)

        # Qual 2 is the default for clair3
        #  CL arg gives options for what we want to keep as a min quality
        if qual < self.min_variant_qual:
            return False

        # Non-divisible by 3 indels are more tolerated at different positions and in different viruses
        #  So allow adjustable min non-divisible qual
        if not in_frame(v):
            if self.no_frameshifts:
                return False
            # Require a higher quality for frameshifting indels, they're far more likely to be errors
            if qual < self.min_frameshift_qual:
                return False

        # Allele frequency
        if allele_freq < self.min_allele_freq:
            return False

        # Depth
        try:
            depth = v.INFO["DP"]
        except KeyError:
            depth = v.format("DP")[0][0]

        if depth < self.min_depth:
            return False

        return True


def go(args):
    vcf_reader = VCF(args.inputvcf)
    vcf_writer = Writer(args.output_pass_vcf, vcf_reader, "w")
    vcf_writer.write_header()
    vcf_writer_filtered = Writer(args.output_fail_vcf, vcf_reader, "w")
    vcf_writer_filtered.write_header()

    if args.nanopolish:
        filter = NanopolishFilter(args.no_frameshifts)
    elif args.medaka:
        filter = MedakaFilter(args.no_frameshifts)
    elif args.clair3:
        filter = Clair3Filter(
            args.no_frameshifts, args.min_depth,
            args.min_qual_c3, args.min_frameshift_qual,
            args.min_allele_freq
        )
    else:
        print("Please specify a VCF type, i.e. --nanopolish or --medaka or --clair3\n")
        raise SystemExit

    variants = [v for v in vcf_reader]

    group_variants = defaultdict(list)
    for v in variants:
        indx = f"{v.CHROM}-{v.POS}"
        group_variants[indx].append(v)

    for v in variants:

        # Pre-filter to remove rubbish that we don't want adding to the mask
        try:
            if v.INFO["DP"] <= 1:
                print(f"Suppress variant {v.POS} due to low depth")
                continue
        except KeyError:
            pass

        # Completely skip RefCalls in clair3
        if v.ALT == []:
            print(f"skipping RefCall at {v.POS}")
            continue

        # No longer focused on medaka but keep the qual prefilter in
        if args.medaka:
            if v.QUAL < 20:
                continue

        # Clair3 prefilter for lowqual and low AF reads to not mask them
        if args.clair3:
            if v.QUAL < 2:
                print(f"Skipping LowQual of {v.QUAL} at {v.POS}")
                continue
            # Skip really low AF to help not get a lot of Ns with noisy data
            allele_freq = v.format("AF")[0][0]
            if allele_freq < args.min_mask_freq:
                print(f"Skipping LowAF of {allele_freq} at {v.POS}")
                continue

        # Now apply the filter to send variants to PASS or FAIL file
        if filter.check_filter(v):
            vcf_writer.write_record(v)
        else:
            variant_passes = False

            indx = f"{v.CHROM}-{v.POS}"
            if len(group_variants[indx]) > 1:
                for check_variant in group_variants[indx]:
                    if filter.check_filter(check_variant):
                        variant_passes = True

            if not variant_passes:
                vcf_writer_filtered.write_record(v)
            else:
                print (f"Suppress variant {v.POS}\n")

def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--nanopolish', action='store_true')
    parser.add_argument('--medaka', action='store_true')
    parser.add_argument('--clair3', action='store_true')
    parser.add_argument('--no-frameshifts', action='store_true')
    parser.add_argument("--min-depth", type=int, default=20)
    parser.add_argument('--min-qual-c3', type=int, default=7)
    parser.add_argument('--min-frameshift-qual', type=int, default=15)
    parser.add_argument('--min-allele-freq', type=float, default=0.60)
    parser.add_argument('--min-mask-freq', type=float, default=0.25)
    parser.add_argument('inputvcf')
    parser.add_argument('output_pass_vcf')
    parser.add_argument('output_fail_vcf')

    args = parser.parse_args()

    go(args)

if __name__ == "__main__":
    main()
