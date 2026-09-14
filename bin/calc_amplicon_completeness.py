#!/usr/bin/env python3
'''Check how complete an amplicon is based bed file and consensus'''
import argparse
import csv
from Bio import SeqIO

def init_parser() -> argparse.ArgumentParser:
    """
    Purpose
    -------
    Parse CL inputs to be used in script

    Returns
    -------
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
        '-a',
        '--amplicon_bed',
        required=True,
        type=str,
        help='Path to amplicon bed file'
    )
    parser.add_argument(
        '-c',
        '--consensus',
        required=True,
        type=str,
        help='Input sample consensus sequence file'
    )
    return parser

def main() -> None:
    '''Run the program'''
    # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    # Read consensus file
    consensus = SeqIO.read(args.consensus, "fasta")
    consensus.seq = consensus.seq.upper()

    # Read in the amplicon bed file
    #  It shouldn't change format as it comes from the primers_to_amplicons script
    with open(args.amplicon_bed) as handle:
        reader = csv.reader(handle, delimiter='\t')
        out = {}
        for row in reader:
            # Remember that by primer conventions start is 0-based and stop is 1-based
            #  So works in python indexes with no changes needed
            start, stop, name = int(row[1]), int(row[2]), str(row[3])
            amp_length = len(range(start, stop))
            n_count = consensus.seq[start:stop].count('N')

            # Calc
            if amp_length == 0:
                continue
            else:
                completeness = round(1 - (n_count/amp_length), 2)
            out[name] = str(completeness)

    # Output
    with open(f'{str(args.sample)}_amplicon_completeness.csv', 'w') as f:
        header = f'sample,{",".join(out.keys())}'
        line = f'{str(args.sample)},{",".join(out.values())}'
        f.write(header)
        f.write('\n')
        f.write(line)
        f.write('\n')

if __name__ == '__main__':
    main()
