#!/usr/bin/env python3
'''Create QC CSV file based on pipeline outputs'''
import argparse
import csv
import statistics
import subprocess
import pandas as pd
import vcf
import re

from Bio import SeqIO, SeqRecord
from collections import defaultdict
from typing import Tuple, Optional

# Global REGEXES
NEXTCLADE_STOP_PATTERN = re.compile(":\\*")


def init_parser() -> argparse.ArgumentParser:
    """Parse CL inputs to be used in script

    Returns:
    --------
        argparse.ArgumentParser
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
        '-i',
        '--irida_id',
        required=True,
        type=str,
        help='IRIDA ID to upload metadata to'
    )
    parser.add_argument(
        '-b',
        '--bam',
        required=True,
        type=str,
        help='Bam file used to call variants from'
    )
    parser.add_argument(
        '-v',
        '--vcf',
        required=True,
        type=str,
        help='Input sample passing vcf file'
    )
    parser.add_argument(
        '-c',
        '--consensus',
        required=True,
        type=str,
        help='Final sample consensus sequence file'
    )
    parser.add_argument(
        '-d',
        '--depth',
        required=True,
        type=str,
        help='Positional depth bed file from samtools depth'
    )
    parser.add_argument(
        '--min_vcf',
        required=False,
        type=str,
        help='Input sample minor vcf file'
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
    parser.add_argument(
        '--add_nextclade_columns',
        action='store_true',
        help='Add nextclade columns to the final output only if params.skip_nextclade is not invoked'
    )
    parser.add_argument(
        '--nextclade_csv',
        required=False,
        type=str,
        help='Nextclade CSV file'
    )
    return parser


def validate_df_columns(df: pd.DataFrame, needed_columns: list) -> None:
    """Check that input CSV contains the correct columns needed. Exits program if not

    Params:
    -------
        df (DataFrame): Dataframe made from the input CSV file
        needed_columns (list): List of columns to confirm exist
    """
    columns = list(df.columns)
    if any(x not in columns for x in needed_columns):
        missing_str = ', '.join([x for x in needed_columns if x not in columns])
        raise ValueError(f'Missing {missing_str} column(s) needed for validation'.format())


def get_read_count(bam: str, chrom: Optional[str] = None) -> int:
    """Get the number of aligned reads from given bamfile and optionally the given chrom using samtools view and subprocess

    Params:
    -------
        bam (str): Path to the bam file used to call variants
        chrom (None | str): Optional name of the segment/chromosome to get the stats for

    Returns:
    --------
        integer: Number of aligned reads
    """
    cmd = ['samtools', 'view', '-c', '-F0x900', bam]
    if chrom:
        cmd.append(chrom)
    read_count = subprocess.run(cmd, capture_output=True, check=True, text=True).stdout.strip('\n')
    return int(read_count)


def parse_depth_bed(bed: str) -> defaultdict:
    """Calculate the mean and median sequencing depth for each chrom from samtools depth bed file

    Params:
    -------
        bed (str): Path to the input depth bed file from samtools depth

    Returns:
    --------
        defaultdict: Containing chrom as keys with median and mode underlying
    """
    depth = defaultdict(list)
    with open(bed) as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for d in reader:
            depth[d['chrom']].append(int(d['depth']))

    # Empty file return nothing
    depth_stats = defaultdict(dict)
    for chrom, depth_list in depth.items():
        if depth_list == []:
            depth_stats[chrom]['mean'] = 0
            depth_stats[chrom]['median'] = 0
            continue
        # Otherwise calc
        mean_dep = round(statistics.mean(depth_list), 2)
        median_dep = round(statistics.median(depth_list), 1)
        depth_stats[chrom]['mean'] = mean_dep
        depth_stats[chrom]['median'] = median_dep
    return depth_stats


def parse_consensus(fasta: SeqRecord) -> Tuple[int, float]:
    """Parse consensus file to get the genome completeness and N count

    Params:
    -------
        fasta (SeqRecord): Single consensus sequence input

    Returns:
    --------
        integer: Number of Ns
        float: Genome completeness calculated from the number of Ns and the genome length
    """
    n_pos =  [i for i, base in enumerate(fasta.seq.lower()) if base == 'n']
    count_n = len(n_pos)
    completeness = 1 - (count_n / len(fasta.seq))
    completeness = round(completeness, 4)
    return count_n, completeness


def _create_variantpos_dict(var: str, var_range: range) -> dict:
    """Create variant position dict for tracking variants

    Params:
    -------
        var (str): Variant
        var_range (range): Genomic range the variant spans

    Returns:
    --------
        dict: with keys "variant" and "range"
    """
    return {
        'variant': var,
        'range': var_range
    }


def parse_vcf(vcf_file: str, chrom: str) -> Tuple[str, list, str, dict]:
    """Parse input VCF file to find variants and their locations

    Params:
    -------
        vcf_file (str): Path to input gzipped vcf file
        chrom (str): Name of chrom to select variants from

    Returns:
    --------
        tuple: Elements parsed from the vcf file for tracking
            str: All parsed variants to report formatted as 'RefPosAlt'
            list: Containing variant-range dicts
            str: Potential frameshift variants (%3==0 check and snpeff ann check)
            dict: Tracking the different variant counts
                {'total_variants': int, 'num_snps': int, 'num_deletions': int, 'num_deletion_sites': int, 'num_insertions': int, 'num_insertion_sites': int}
    """
    # Base outputs and counts
    variants = []
    variant_positions = []
    frameshift_variants = []
    var_count_dict = {
        'total_variants': 0,
        'num_snps': 0,
        'num_iupacs': 0,
        'num_deletions': 0,
        'num_deletion_sites': 0,
        'num_insertions': 0,
        'num_insertion_sites': 0
    }

    # Open and handle file
    with open(vcf_file, 'rb') as handle:
        reader = vcf.Reader(handle)
        for record in reader:
            # Only wanted chrom allowed
            if record.CHROM != chrom:
                continue

            # Base information
            ref_len = len(record.REF)
            alt_len = len(record.ALT[0])
            alt_str = str(record.ALT[0])
            iupac = False
            if ("ConsensusTag" in record.INFO) and (record.INFO["ConsensusTag"] == 'ambiguous'):
                iupac = True
                alt_str = record.INFO["ConsensusBase"]

            # Multiple alleles should be removed by now by all methods. Exit if not for bugfixing
            #  IUPAC do contain multiple though so that is the only exception here
            if len(record.ALT) > 1:
                if not iupac:
                    raise ValueError(f"Multiple alleles should have been resolved previously. Check intermediate VCFs at position: {record.POS}")
            # Odd issue previously - skip over Ns in vcf
            if str(record.ALT[0]).upper() == 'N':
                continue
            # Other odd issue, skip positions where alt is None
            if record.ALT[0] is None:
                continue

            # Create string of variant and add to list of variants along with getting the lengths of the REF and ALT
            variant = f'{record.REF}{record.POS}{alt_str}'

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
                    # Again, the alt will include 1 reference position so minus 1
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
                #  Note: Should have been resolved previously
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
                    if iupac:
                        var_count_dict['num_iupacs'] += 1
                    else:
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


def range_contains(r1: range, r2: range) -> bool:
    """Returns True if range1 fully contains range2"""
    x1, x2 = r1.start, r1.stop
    y1, y2 = r2.start, r2.stop
    return x1 <= y2 and y1 <= x2


def check_primers(bed: str, variant_locations: list, chrom: str) -> str:
    """Parse input bed file for any variant regions that overlap

    Params:
    -------
        bed (str): Path to primer bed file
        variant_locations (list[dict]): Variant position dictionaries with variant and range keys
        chrom (str): Name of the segment the variants are from

    Returns:
    --------
        str: Primer mutations found or string "none" if there were none
    """
    if not variant_locations:
        return 'none'

    primer_mutations = []
    with open(bed) as handle:
        reader = csv.reader(handle, delimiter='\t')
        for row in reader:
            # Bed file needs at least 4 rows (chrom, start, stop, name)
            if len(row) < 4:
                continue
            # Only check on the right chrom for variants
            if row[0] != chrom:
                continue

            # Set primer values, make sure start lower than stop for range
            start, stop, name = int(row[1]), int(row[2]), str(row[3])
            if start > stop:
                start, stop = stop, start
            location = range(start, stop + 1) # Plus one to make sure that we get mutations in the final location of the range

            # Check if the location range overlaps with any variant ranges
            for var_dict in variant_locations:
                if range_contains(location, var_dict['range']):
                    primer_mutations.append(f'{var_dict["variant"]}-{name}')

    if primer_mutations:
        return ';'.join(primer_mutations)
    return 'none'


def parse_metadata(metadata: str, sample: str) -> pd.DataFrame:
    """Parse metadata file to find metadata for given sample

    Params:
    -------
        metadata (str): Path to metadata file
        sample  (str): Sample name to look for

    Returns:
    --------
        DataFrame: Containing columns from the wanted sample or empty df
    """
    df = pd.read_csv(metadata, sep='\t')
    validate_df_columns(df, ['sample'])
    df = df.loc[df['sample'] == sample]
    if len(df) == 1:
        return df
    elif len(df) > 1:
        raise RuntimeError(f'Sample {sample} exists more than once in the metadata file')
    else:
        # Can return an empty df
        return df


def count_minor_variants(vcf_file: str, chrom: str) -> Tuple[int, int]:
    """Small function to count passing SNPs and indels in the minor VCF file.

    Params:
    -------
        vcf_file (str): Path to the minor VCF file.
        chrom (str): Chromosome to filter variants.

    Returns:
    --------
        Tuple[int, int]: Number of passing SNPs and indels.
    """
    snps = 0
    indels = 0

    with open(vcf_file, 'rb') as handle:
        reader = vcf.Reader(handle)
        for record in reader:
            if record.CHROM != chrom:
                continue
            if record.FILTER and "PASS" not in record.FILTER:
                continue
            ref_len = len(record.REF)
            alt_len = len(record.ALT[0])
            if ref_len == 1 and alt_len == 1:
                snps += 1
            elif ref_len != alt_len:
                indels += 1

    return snps, indels


def get_nextclade_vals(nextclade_csv: str) -> Tuple[str, str, str, int]:
    '''Parse custom nextclade CSV file to find information on potential issue sites

    Parameters:
    -----------
        nextclade_csv (str): Path to nextclade CSV file. ';' delimited

    Returns:
    --------
        Tuple[str, str, str, int] : frameshifts, stop codons, mutated stop codons, frameshift count
    '''
    # Nextclade CDS Checks
    ## To look back at for segmented viruses since it now just uses the second line
    with open(nextclade_csv) as handle:
        reader = csv.DictReader(handle, delimiter=';')
        d = next(reader, None)

    if d:
        aa_mutations = d['aaSubstitutions']
        frameshifts = d['qc.frameShifts.frameShifts']
        stop_codons = d['qc.stopCodons.stopCodons']
        # Failing samples or unmatched samples have the column but as an empty string so need the fallback 0
        total_fs = d['qc.frameShifts.totalFrameShifts'] or 0
        ignored_fs = d['qc.frameShifts.totalFrameShiftsIgnored'] or 0

        mutated_stop_codons_match = re.findall(NEXTCLADE_STOP_PATTERN, aa_mutations)
        mutated_stop_codons = '|'.join(mutated_stop_codons_match)

        # Have to set these to int as dict reader is bringing the values back as strings
        fs_count = int(total_fs) - int(ignored_fs)

        return frameshifts, stop_codons, mutated_stop_codons, fs_count
    else:
        return '', '', '', 0


def grade_qc(completeness: float, mean_dep: float, median_dep: float,
             fs_status: bool, nonsense_status: bool, mutated_stop_status: bool) -> str:
    """Determine if the sample passes internal QC metrics and assign a PASS or why it failed

    Params:
    -------
        completeness (float): Final genome completeness
        mean_dep (float): Mean sequencing depth
        median_dep (float): Median sequencing depth
        fs_status (bool): True if any potential frameshift variants from VCF parsing or Nextclade are found
        nonsense_status (bool): True if any nonsensus mutations from Nextclade are found
        mutated_stop_status (bool): True if any stop codon mutations from Nextclade are found

    Returns:
    --------
        str: Final QC status
    """
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
    if fs_status:
        qc_status.append('POTENTIAL_FRAMESHIFTS')
    # Nonsense
    if nonsense_status:
        qc_status.append('NONSENSE_MUTATION')
    # Stop Codon Mutations
    if mutated_stop_status:
        qc_status.append('STOP_CODON_MUTATION')

    if qc_status:
        return ';'.join(qc_status)
    return 'PASS'


def main() -> None:
    """Main entry to the program"""
    # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    # Checks that Chrom in them already and output a dict or overall setting
    #  Depth/Reads
    depth_dict = parse_depth_bed(args.depth)
    total_reads = get_read_count(args.bam)

    # Per chrom/segment information and final output setup
    final_out = []
    with open(args.consensus) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # The chrom/segment is ALWAYS after the sample name based on how the pipeline is setup to rename consensus seqs
            #  and there are no spaces allowed in the sample name
            chrom = record.description.split(' ')[1]

            # Reads
            num_reads = get_read_count(args.bam, chrom)
            count_n, completeness = parse_consensus(record)

            # Variants
            variants, variant_positions, frameshift_variants, var_count_dict = parse_vcf(args.vcf, chrom)

            # Optional primer checks, same as variants parsing multiple times for now
            pcr_primer_overlap = 'NA'
            seq_primer_overlap = 'NA'
            if args.seq_bed:
                seq_primer_overlap = check_primers(args.seq_bed, variant_positions, chrom)
            if args.pcr_bed:
                pcr_primer_overlap = check_primers(args.pcr_bed, variant_positions, chrom)

            # Minor variants (if provided)
            if args.min_vcf:
                minor_snps, minor_indels = count_minor_variants(args.min_vcf, chrom)

            # Nextclade mutations (if provided)
            nc_frameshifts, nc_nonsense, nc_mutated_stop_codons, nc_fs_count = '', '', '', 0
            if args.add_nextclade_columns and args.nextclade_csv:
               nc_frameshifts, nc_nonsense, nc_mutated_stop_codons, nc_fs_count = get_nextclade_vals(args.nextclade_csv)

            # Grade QC
            mean_depth = depth_dict[chrom].get('mean', 0)
            median_depth = depth_dict[chrom].get('median', 0)
            #  Using the nextclade count (if available) to take into account the full gene effect
            fs_status = ((frameshift_variants != 'none') and (nc_fs_count > 0))
            nonsense_status = (nc_nonsense != '')
            mutated_stop_status = (nc_mutated_stop_codons != '')

            qc_status = grade_qc(completeness, mean_depth, median_depth, fs_status, nonsense_status, mutated_stop_status)

            # Final Output
            #  Main columns
            sample_data = {
                'sample': args.sample,
                'reference': chrom,
                'qc_pass': qc_status,
                'num_aligned_reads': total_reads,
                'num_segment_reads': num_reads,
                'num_consensus_n': count_n,
                'genome_completeness': completeness,
                'mean_sequencing_depth': mean_depth,
                'median_sequencing_depth': median_depth,
                'total_variants': var_count_dict['total_variants'],
                'num_snps': var_count_dict['num_snps'],
                'num_iupacs': var_count_dict['num_iupacs'],
                'num_deletions': var_count_dict['num_deletions'],
                'num_deletion_sites': var_count_dict['num_deletion_sites'],
                'num_insertions': var_count_dict['num_insertions'],
                'num_insertion_sites': var_count_dict['num_insertion_sites']
            }

            # Conditionally add nextclade mutation data
            if args.add_nextclade_columns:
                sample_data['nextclade_frameshifts'] = nc_frameshifts
                sample_data['nonsense_mutations'] = nc_nonsense
                sample_data['mutated_stop_codons'] = nc_mutated_stop_codons

            # Conditionally add the minor variant data
            if args.min_vcf:
                sample_data['minor_snps'] = minor_snps
                sample_data['minor_indels'] = minor_indels

            # Remaining 'busy' columns
            sample_data['variants'] = variants
            sample_data['possible_frameshift_variants'] = frameshift_variants
            sample_data['sequencing_primer_variants'] = seq_primer_overlap
            sample_data['diagnostic_primer_variants'] = pcr_primer_overlap
            sample_data['irida_id'] = args.irida_id

            final_out.append(sample_data)


    # Create and output final CSV
    df = pd.DataFrame.from_dict(final_out)
    if args.metadata:
        metadata_df = parse_metadata(args.metadata, args.sample)
        df = df.merge(metadata_df, on='sample', how='left')
    df.to_csv(f'{args.sample}.qc.csv', index=False)


if __name__ == '__main__':
    main()
