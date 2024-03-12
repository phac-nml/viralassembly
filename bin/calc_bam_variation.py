#!/usr/bin/env python3

import argparse
import pysam
import pandas as pd
import sys
from pathlib import Path

def init_parser() -> argparse.ArgumentParser:
    '''
    Purpose
    -------
    Parse CL inputs to be used in script
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-b',
        '--bam',
        required=True,
        help='Path to input BAM file. Make sure that the BAI file is also on the same path'
    )
    parser.add_argument(
        '-r',
        '--reference',
        required=True,
        help='Input reference fasta file'
    )
    parser.add_argument(
        '--sample',
        required=False,
        help='Input sample name to be used in output file creation'
    )
    parser.add_argument(
        '--min_base_qual',
        required=False,
        type=int,
        default=10,
        help='Minimum base quality to use read for position'
    )
    parser.add_argument(
        '--min_map_qual',
        required=False,
        type=int,
        default=30,
        help='Minimum mapping quality to use read for position'
    )
    parser.add_argument(
        '--min_read_count',
        required=False,
        type=int,
        default=20,
        help='Minimum number of reads to report position'
    )
    parser.add_argument(
        '--min_report_percent',
        required=False,
        type=int,
        default=15,
        help='Minimum percentage of non-reference bases to report out'
    )
    parser.add_argument(
        '--variant_call_threshold',
        required=False,
        type=int,
        default=70,
        help='Minimum percentage for a variant to called'
    )

    return parser

def create_ref_dict(ref_fasta) -> dict:
    """
    Purpose
    -------
    Process reference fasta file with pysam to create a dictionary of indexs with their bases

    Parameters
    ----------
    ref_fasta : str
        Path to input reference fasta file

    Returns
    -------
    dict{ Position: Base }
    """
    ref_dict = {}
    ref_fasta = pysam.FastaFile(ref_fasta)
    for idx, base in enumerate(ref_fasta.fetch(ref_fasta.references[0])):
        ref_dict[idx] = base

    return ref_dict

def calculate_non_ref_percent(pos, pileup_dict, ref_dict) -> float:
    """
    Purpose
    -------
    Calculate the non-reference base percentage for the position

    Parameters
    ----------
    pos : int
        Current genome position
    pileup_dict: dict
        Contains read counts for the current position
    ref_dict: dict
        Contains: {Position: Base} for the reference genome

    Returns
    -------
    Non-reference percentage as float
    """
    total = pileup_dict["total_reads"]
    if ref_dict[pos] != "N":
        ref_base = ref_dict[pos] + "_reads"
        ref_count = pileup_dict[ref_base]
        if (total == 0):
            non_ref_percent = 0
        elif (ref_count == 0):
            non_ref_percent = 100
        else:
            non_ref_percent = round((100 - ((ref_count / total) * 100)), 2)
    else:
        non_ref_percent = 0

    return non_ref_percent

def calc_base_percent(pileup_dict, base) -> float:
    """
    Purpose
    -------
    Calculate the base percentage for the position

    Parameters
    ----------
    pileup_dict: dict
        Contains read counts for the current position
    base: str
        Dict key of the base to use for percent calculation

    Returns
    -------
    Given base percentage as float
    """
    if pileup_dict["total_reads"] == 0:
        return 0
    return round(((pileup_dict[base] / pileup_dict["total_reads"]) * 100), 2)

def which_majority(pileup_dict, variant_call_threshold) -> bool:
    """
    Purpose
    -------
    Determine what the called mutation type is for the given pileup_dict

    Parameters
    ----------
    pileup_dict: dict
        Contains read counts for the current position
    variant_call_threshold: int
        Threshold at which a variant is called

    Returns
    -------
    String of mutation type
    """
    if pileup_dict['A_percent'] >= variant_call_threshold:
        return 'A SNP'
    elif pileup_dict['C_percent'] >= variant_call_threshold:
        return 'C SNP'
    elif pileup_dict['T_percent'] >= variant_call_threshold:
        return 'T SNP'
    elif pileup_dict['G_percent'] >= variant_call_threshold:
        return 'G SNP'
    else:
        return 'Del'

def determine_type(pileup_dict, variant_call_threshold) -> str:
    """
    Purpose
    -------
    Determine if the location is a mutation, mixed site, or reference base

    Parameters
    ----------
    ref_dict: dict
        Contains: {Position: Base} for the reference genome
    variant_call_threshold: int
        Threshold at which a variant is called. Default: 70

    Returns
    -------
    String of the location status
    """
    if pileup_dict['percentage_nonref'] >= variant_call_threshold:
        return which_majority(pileup_dict, variant_call_threshold)
    elif (100 - variant_call_threshold) <= pileup_dict['percentage_nonref'] < variant_call_threshold:
        return 'Mixed'
    else:
        return 'Ref'

def parse_variation_from_bam(bamfile, ref_dict, base_q, map_q, min_read_count, min_report_percent, variant_call_threshold) -> list:
    """
    Purpose
    -------
    Parse BAM file using pysam mpileup to find positions where non-ref base count is above
        the min_report_percent percent and min_read_count

    Parameters
    ----------
    bamfile : str
        Path to input BAM file
    ref_dict: dict
        Contains: {Position: Base} for the reference genome
    base_q: int
        Minimum base quality to use read. Default: 10
    map_q: int
        Minimum mapping quality to use read. Default: 30
    min_read_count: int
        Minimum reads required to report position. Default: 20
    min_report_percent: int
        Minimum non-reference base percentage to report out. Default: 15
    variant_call_threshold: int
        Threshold at which a variant is called. Default: 70

    Returns
    -------
    list of dictionaries containing positional variation information if variation is above min_report_percent
    """
    bamfile = pysam.AlignmentFile(bamfile, "rb" )
    variation_info = []

    # Based on pysam docs for iterating over each base and Piranha's variation analysis
    for pileupcolumn in bamfile.pileup(bamfile.references[0], min_base_quality=base_q, min_mapping_quality=map_q):
        pileup_dict = {}
        A_counter, G_counter, C_counter, T_counter, del_counter = 0,0,0,0,0

        # For each read in the specific position, count total base found
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                # Query position is None if is_del or is_refskip is set
                if pileupread.alignment.query_sequence[pileupread.query_position] == "A":
                    A_counter += 1
                elif pileupread.alignment.query_sequence[pileupread.query_position] == "G":
                    G_counter += 1
                elif pileupread.alignment.query_sequence[pileupread.query_position] == "C":
                    C_counter += 1
                elif pileupread.alignment.query_sequence[pileupread.query_position] == "T":
                    T_counter += 1
            else:
                del_counter += 1

        # Add values to dictionary
        #  chrom is second as multiqc uses first column for output and having it be the same leads to issues
        pileup_dict["position"] = pileupcolumn.reference_pos + 1
        pileup_dict["chrom"] = pileupcolumn.reference_name
        pileup_dict["A_reads"] = A_counter
        pileup_dict["C_reads"] = C_counter
        pileup_dict["T_reads"] = T_counter
        pileup_dict["G_reads"] = G_counter
        pileup_dict["del_reads"] = del_counter
        pileup_dict["total_reads"] = (A_counter + C_counter + T_counter + G_counter + del_counter)
        pileup_dict["percentage_nonref"] = calculate_non_ref_percent(pileupcolumn.pos, pileup_dict, ref_dict)
        pileup_dict["ref_base"] = ref_dict[pileupcolumn.reference_pos]
        pileup_dict["A_percent"] = calc_base_percent(pileup_dict, "A_reads")
        pileup_dict["C_percent"] = calc_base_percent(pileup_dict, "C_reads")
        pileup_dict["T_percent"] = calc_base_percent(pileup_dict, "T_reads")
        pileup_dict["G_percent"] = calc_base_percent(pileup_dict, "G_reads")
        pileup_dict["del_percent"] = calc_base_percent(pileup_dict, "del_reads")
        pileup_dict["variant_type"] = determine_type(pileup_dict, variant_call_threshold)

        # Report based on threshold and deletion status
        if pileup_dict["percentage_nonref"] >= min_report_percent:
            if pileup_dict["total_reads"] >= min_read_count:
                variation_info.append(pileup_dict)

    return variation_info

def main() -> None:
    '''
    Main entry point into script
    '''
    ## Init Parser and Set Arguments ##
    parser = init_parser()
    args = parser.parse_args()

    # Current test analysis
    path_bam = Path(args.bam)
    if args.sample:
        sample_name = args.sample
    else:
        sample_name = path_bam.name.split('.')[0]
    if path_bam.stat().st_size < 200:
        print(f"{args.bam} is too small a file. Skipping")
        sys.exit(0)

    ref_dict = create_ref_dict(args.reference)
    var_info = parse_variation_from_bam(path_bam, ref_dict, args.min_base_qual,
                                        args.min_map_qual, args.min_read_count,
                                        args.min_report_percent, args.variant_call_threshold)
    if var_info == []:
        sys.exit(0)
    df = pd.DataFrame.from_dict(var_info)
    df.to_csv(f'{sample_name}_variation.csv', index=False)

if __name__ == "__main__":
    main()
