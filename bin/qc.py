#!/usr/bin/env python3
'''Create QC CSV file based on pipeline outputs'''
import argparse
import csv
import statistics
import subprocess
import pandas as pd
import vcf

from Bio import SeqIO, SeqRecord
from typing import Tuple

def init_parser() -> argparse.ArgumentParser:
    """
    Specify command line arguments
    Returns command line parser with inputs
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-s',
        '--sample',
        required=True,
        type=str,
        help='Sample name'
    )
    parser.add_argument(
        '-a',
        '--analysis',
        required=True,
        type=str,
        help='Analysis program used ["nanopolish", "medaka", "clair3"]'
    )
    parser.add_argument(
        '-c',
        '--consensus',
        required=True,
        type=str,
        help='Input sample consensus sequence file'
    )
    parser.add_argument(
        '-b',
        '--bam',
        required=True,
        type=str,
        help='Input sample bam file'
    )
    parser.add_argument(
        '-d',
        '--depth',
        required=True,
        type=str,
        help='Input sample depth bed file'
    )
    parser.add_argument(
        '-v',
        '--vcf',
        required=True,
        type=str,
        help='Input sample passing vcf file'
    )
    parser.add_argument(
        '-m',
        '--metadata',
        required=False,
        type=str,
        help='Input run TSV metadata file'
    )
    parser.add_argument(
        '--seq_bed',
        required=False,
        type=str,
        help='Input Sequencing primer bed file to test variants against'
    )
    parser.add_argument(
        '--pcr_bed',
        required=False,
        type=str,
        help='Input PCR bed file to test variants against'
    )
    return parser

def validate_df_columns(df: pd.DataFrame, needed_columns: list) -> None:
    """
    Purpose
    -------
    Check that input CSV contains the correct columns needed. Exits program if not

    Parameters
    ----------
    df: pd.DataFrame
        Pandas dataframe made from the input CSV file
    """
    columns = list(df.columns)
    if any(x not in columns for x in needed_columns):
        raise ValueError('Missing {} column(s) needed for validation'.format([x for x in needed_columns if x not in columns]))

def get_read_count(bam: str) -> int:
    '''
    Purpose:
    --------
    Get the number of aligned reads from given bamfile using samtools and subprocess

    Parameters:
    -----------
    bam - str
        Path to the primertrimmed bam file

    Returns:
    --------
    Integer number of aligned reads
    '''
    cmd = ['samtools', 'view', '-c', '-F0x900', bam]
    read_count = subprocess.run(cmd, capture_output=True, check=True, text=True).stdout.strip('\n')
    return int(read_count)

def parse_depth_bed(bed: str) -> Tuple[float, float]:
    '''
    Purpose:
    --------
    Calculate the mean and median sequencing depth from samtools depth bed file

    Parameters:
    -----------
    bed - str
        Path to the input depth bed file from samtools depth

    Returns:
    --------
    Float mean sequencing depth
    Float median sequencing depth
    '''
    depth = []
    with open(bed, 'r') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for d in reader:
            depth.append(int(d['depth']))

    # Empty file return nothing
    if depth == []:
        return 0, 0
    # Otherwise calc
    mean_dep = round(statistics.mean(depth), 2)
    median_dep = round(statistics.median(depth), 1)
    return mean_dep, median_dep

def parse_consensus(fasta: SeqRecord) -> Tuple[int, float]:
    '''
    Purpose:
    --------
    Parse consensus file to get the genome completeness and N count

    Parameters:
    -----------
    fasta - SeqRecord
        Bio SeqRecord of input consensus sequence

    Returns:
    --------
    Integer number of Ns
    Float genome completeness calculated from the number of Ns
    '''
    n_pos =  [i for i, base in enumerate(fasta.seq.lower()) if base == 'n']
    count_n = len(n_pos)
    completeness = 1 - (count_n / len(fasta.seq))
    completeness = round(completeness, 4)
    return count_n, completeness

def _create_variantpos_dict(var: str, var_range: range) -> dict:
    '''Create dict with keys "variant" and "range"'''
    return {
        'variant': var,
        'range': var_range
    }

def parse_vcf(vcf_file: str) -> Tuple[str, list, str, str, dict]:
    '''
    Purpose:
    --------
    Parse input VCF file to find variants and their locations

    Parameters:
    -----------
    vcf_file - str
        Path to input gzipped vcf file from args

    Returns:
    --------
    String of parsed variants to report
    List of variant-range dicts
    String of any potential frameshift variants (%3==0 check and snpeff ann check)
    Dict of different variant counts
        {'total_variants': int, 'num_snps': int, 'num_deletions': int, 'num_deletion_sites': int, 'num_insertions': int, 'num_insertion_sites': int}
    '''
    # Base outputs
    variants = []
    variant_positions = []
    frameshift_variants = []
    var_count_dict = {
        'total_variants': 0,
        'num_snps': 0,
        'num_deletions': 0,
        'num_deletion_sites': 0,
        'num_insertions': 0,
        'num_insertion_sites': 0
    }

    # Open and handle file
    with open(vcf_file, 'rb') as handle:
        reader = vcf.Reader(handle)
        for record in reader:
            # Odd issue previously - skip over Ns in vcf
            if str(record.ALT[0]).upper() == 'N':
                continue
            # Other odd issue, skip positions where alt is None
            if record.ALT[0] == None:
                continue
            # Multiple alleles not supported warning
            if len(record.ALT) > 1:
                print(f'WARNING: Multiple alleles not supported currently. Using only the first one for position: {record.POS}')

            # Create string of variant and add to list of variants along with getting the lengths of the ref and alt
            variant = f'{record.REF}{record.POS}{record.ALT[0]}'
            ref_len = len(record.REF)
            alt_len = len(record.ALT[0])
            alt_str = str(record.ALT[0])

            # Type of mutation leads to different spots affected and tracking
            if variant not in variants:
                # Deletions
                if ref_len > alt_len:
                    # For dels, POS is the kept genomic position so +1 to start to get the deleted positions
                    # The range works as a 9bp deletion would be length of 10 in vcf record
                    del_range = range(record.POS+1,record.POS+ref_len)
                    if (len(del_range) %3 != 0):
                        # If we have annotations check those as well
                        var_ann = record.INFO.get('ANN', '')
                        if var_ann:
                            ann_list = var_ann[0].split('|')
                            consequence = str(ann_list[1]).lower()
                            if 'frameshift' in consequence:
                                frameshift_variants.append(variant)
                        else:
                            frameshift_variants.append(variant)
                    variants.append(variant)
                    variant_positions.append(_create_variantpos_dict(variant, del_range))
                    var_count_dict['num_deletions'] += len(del_range)
                    var_count_dict['num_deletion_sites'] += 1

                # Insertions
                elif ref_len < alt_len:
                    if ((alt_len-1) %3 != 0):
                        # If we have annotations check those as well
                        var_ann = record.INFO.get('ANN', '')
                        if var_ann:
                            ann_list = var_ann[0].split('|')
                            consequence = str(ann_list[1]).lower()
                            if 'frameshift' in consequence:
                                frameshift_variants.append(variant)
                        else:
                            frameshift_variants.append(variant)
                    variants.append(variant)
                    variant_positions.append(_create_variantpos_dict(variant, range(record.POS, record.POS+1)))
                    var_count_dict['num_insertions'] += (alt_len - 1) # -1 for the included ref base
                    var_count_dict['num_insertion_sites'] += 1
                # Multiple SNPs together
                elif (ref_len > 1) and (ref_len == alt_len):
                    mult_snp_range = range(record.POS, record.POS+len(record.REF))
                    for i, ref_base in enumerate(record.REF):
                        # Include only snps incase the variant is like ref=ATG alt=TTC where the T isn't a SNP
                        if ref_base == alt_str[i]:
                            pass
                        variant = f'{ref_base}{mult_snp_range[i]}{alt_str[i]}'
                        variants.append(variant)
                        variant_positions.append(_create_variantpos_dict(variant, range(mult_snp_range[i], mult_snp_range[i]+1)))
                        var_count_dict['num_snps'] += 1
                else:
                    variants.append(variant)
                    variant_positions.append(_create_variantpos_dict(variant, range(record.POS, record.POS+1)))
                    var_count_dict['num_snps'] += 1

    # Final Summary and Return
    if frameshift_variants:
        frameshift_variants = ';'.join(frameshift_variants)
    else:
        frameshift_variants = 'none'
    if variants:
        var_count_dict['total_variants'] = len(variants)
        return ';'.join(variants),  variant_positions, frameshift_variants, var_count_dict
    return 'none', variant_positions, frameshift_variants, var_count_dict

def range_overlap(r1: range, r2: range) -> bool:
    '''Return True if range2 overlaps range1'''
    x1, x2 = r1.start, r1.stop
    y1, y2 = r2.start, r2.stop
    return x1 <= y2 and y1 <= x2

def check_primers(bed, variant_locations: list) -> str:
    '''
    Purpose:
    --------
    Parse input bed file for any variant regions that overlap

    Parameters:
    -----------
    bed - str
        Path to bed file
    variant_locations - list[dict]
        List of variantpos dictionaries with variant and range keys

    Returns:
    --------
    String of primer mutations found or string "none" if there are none
    '''
    if not variant_locations:
        return 'none'

    primer_mutations = []
    with open(bed, 'r') as handle:
        reader = csv.reader(handle, delimiter='\t')
        for row in reader:
            # Bed file needs at least 4 rows (chrom, start, stop, name)
            if len(row) < 4:
                continue
            # Set primer values, make sure start lower than stop for range
            start, stop, name = int(row[1]), int(row[2]), row[3]
            if start > stop:
                start, stop = stop, start
            location = range(start, stop + 1) # Plus one to make sure that we get mutations in the final location of the range

            # Check if the location range overlaps with any variant ranges
            for var_dict in variant_locations:
                if range_overlap(location, var_dict['range']):
                    primer_mutations.append(f'{var_dict["variant"]}-{name}')

    if not primer_mutations:
        return 'none'
    return ';'.join(primer_mutations)

def parse_metadata(metadata: str, sample: str) -> pd.DataFrame:
    '''
    Purpose:
    --------
    Parse metadata file for given sample

    Parameters:
    -----------
    metadata - str
        Path to metadata file
    sample - str
        Sample name to look for

    Returns:
    --------
    Dataframe containing the found sample
    '''
    df = pd.read_csv(metadata, sep='\t')
    validate_df_columns(df, ['sample'])
    df = df.loc[df['sample'] == sample]
    if len(df) == 1:
        return df
    elif len(df) > 1:
        raise RuntimeError(f'Sample {sample} exists more than once in the metadata file')
    else:
        # Can probably return an empty df?
        return df

def grade_qc(completeness: float, mean_dep: float, median_dep: float, frameshift_vars: str) -> str:
    '''
    Purpose:
    --------
    Determine if the sample passes internal QC metrics

    Parameters:
    -----------
    completeness - float
        Genome completeness
    mean_dep - float
        Mean sequencing depth
    median_dep - float
        Median sequencing depth
    frameshift_vars - str
        String of any potential frameshift variants or "none" if there weren't any

    Returns:
    --------
    String QC status
    '''
    qc_status = []
    # Completeness
    if completeness < 0.9:
        if completeness < 0.5:
            qc_status.append('INCOMPLETE_GENOME')
        else:
            qc_status.append('PARTIAL_GENOME')
    # Coverage Depth
    if (mean_dep < 20) or (median_dep < 20):
        qc_status.append('LOW_SEQ_DEPTH')
    # Frameshifts
    if frameshift_vars != 'none':
        qc_status.append('POTENTIAL_FRAMESHIFTS')

    if qc_status:
        return ';'.join(qc_status)
    return 'PASS'

def main() -> None:
    '''Run the program'''
    # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    # Do something
    num_reads = get_read_count(args.bam)
    mean_dep, median_dep = parse_depth_bed(args.depth)
    consensus = SeqIO.read(args.consensus, "fasta")
    count_n, completeness = parse_consensus(consensus)
    variants, variant_positions, frameshift_variants, var_count_dict = parse_vcf(args.vcf)

    # Optional inputs
    if args.metadata:
        metadata_df = parse_metadata(args.metadata, args.sample)
    pcr_primer_overlap = 'NA'
    seq_primer_overlap = 'NA'
    if args.seq_bed:
        seq_primer_overlap = check_primers(args.seq_bed, variant_positions)
    if args.pcr_bed:
        pcr_primer_overlap = check_primers(args.pcr_bed, variant_positions)

    # Grade qc
    qc_status = grade_qc(completeness, mean_dep, median_dep, frameshift_variants)

    # Output
    final = {
        'sample': [args.sample],
        'num_aligned_reads': [num_reads],
        'num_consensus_n': [count_n],
        'genome_completeness': [completeness],
        'mean_sequencing_depth': [mean_dep],
        'median_sequencing_depth': [median_dep],
        'total_variants': [var_count_dict['total_variants']],
        'num_snps': [var_count_dict['num_snps']],
        'num_deletions': [var_count_dict['num_deletions']],
        'num_deletion_sites': [var_count_dict['num_deletion_sites']],
        'num_insertions': [var_count_dict['num_insertions']],
        'num_insertion_sites': [var_count_dict['num_insertion_sites']],
        'variants': [variants],
        'possible_frameshift_variants': [frameshift_variants],
        'sequencing_primer_variants': [seq_primer_overlap],
        'diagnostic_primer_variants': [pcr_primer_overlap],
        'qc_pass': [qc_status]
    }
    df = pd.DataFrame.from_dict(final)
    if args.metadata:
        df = df.merge(metadata_df, on='sample', how='left')
    df.to_csv(f'{args.sample}.qc.csv', index=False)


if __name__ == '__main__':
    main()
