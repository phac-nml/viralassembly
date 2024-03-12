#!/usr/bin/env python3
'''
Simple script to parse snpeff annotated VCF file to create a tsv file
    Based off of ncov-tools convert_recurrent_snpeff_data.py script
'''

import argparse
import vcf

def init_parser() -> argparse.ArgumentParser:
    '''
    Specify command line arguments
    Returns command line parser with inputs
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-s',
        '--sample',
        required=True,
        type=str,
        help='Sample name'
    )
    parser.add_argument(
        '-v',
        '--vcf',
        required=True,
        type=str,
        help='Path to annotated VCF file'
    )
    parser.add_argument(
        '-o',
        '--outfile',
        required=False,
        type=str,
        help='Name of outfile. Default is "SAMPLE.vcf.tsv"'
    )
    parser.add_argument(
        '--annotated',
        required=False,
        default=False,
        action='store_true',
        help='Specify if the vcf is annotated or not'
    )
    return parser

def process_variants_details(var, annotated: bool, variants_analyzed=[]) -> dict:
    """
    Process variant and create dictionary for entry
    Remove duplicated variants
    """
    # Check if seen before - if not add it to list and continue on
    variant_str = f'{var.REF}{var.POS}{var.ALT[0]}' # Only first alt allele, shouldn't have more than one with the process currently
    if variant_str in variants_analyzed:
        return {}
    variants_analyzed.append(variant_str)

    # Non-annotated we just need the info in TSV format
    if not annotated:
        out = {
            'Pos': int(var.POS),
            'Ref': str(var.REF),
            'Alt': str(var.ALT[0]),
            'Qual': str(var.QUAL)
        }
        return out

    # VCF ANN is formatted as a "," list to start so if we find it, we take the first entry and split on |
    #  Also want to make sure that we don't have duplicates by checking if our variant str has been seen already
    var_ann = var.INFO.get('ANN', '')
    if not var_ann:
        out = {
            'Pos': int(var.POS),
            'Variant': '',
            'Consequence': '',
            'Ref': str(var.REF),
            'Alt': str(var.ALT[0]),
            'Qual': str(var.QUAL)
        }
        return out
    ann_list = var_ann[0].split('|')

    # Create output detail dict based on VCFannotation format position
    #  https://pcingola.github.io/SnpEff/adds/VCFannotationformat_v1.0.pdf
    impact = ann_list[2]
    gene = ann_list[3]
    prot = ann_list[10]
    pro_var = f'{gene}'
    if prot:
        pro_var = f'{gene} {prot}'
    consequence = ann_list[1]
    out = {
        'Pos': int(var.POS),
        'Variant': pro_var,
        'Consequence': consequence,
        'Ref': str(var.REF),
        'Alt': str(var.ALT[0]),
        'Qual': str(var.QUAL)
    }
    return out

def write_outfile(outfile: str, variants: list) -> None:
    """Write output file to given location"""
    columns = variants[0].keys()
    header = '\t'.join(columns)
    with open(outfile, 'w') as f:
        f.write(header)
        f.write('\n')
        for var_dict in variants:
            f.write('\t'.join([str(var_dict[col]) for col in columns]))
            f.write('\n')

def main() -> None:
    """Main process"""
    # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    # Load in VCF and parse
    vcf_reader = vcf.Reader(filename=args.vcf)
    variants = []
    # To fix out of range issue from pyvcf on empty files
    try:
        for var in vcf_reader:
            var_dict = process_variants_details(var, args.annotated)
            if var_dict != {}:
                variants.append(var_dict)
    except IndexError:
        pass

    # Write outfile
    outfile = f'{args.sample}.vcf.tsv'
    if args.outfile:
        outfile = args.outfile

    if variants != []:
        write_outfile(outfile, variants)
    else:
        print("No Variants")

if __name__ == '__main__':
    main()
