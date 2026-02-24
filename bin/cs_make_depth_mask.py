#!/usr/bin/env python3
"""
Script to create a genomic depth mask for consensus generation
    Adapted from https://github.com/artic-network/fieldbioinformatics/blob/master/artic/make_depth_mask.py
    to remove the need for RG tags steps for non-amplicon data and some small reformats
"""
from Bio import SeqIO
import itertools
import os
import pysam


def collect_depths(bamfile, ref_name, min_depth, ignore_deletions):
    """Collect read depth of coverage per reference position in a BAM file.

    Parameters
    ----------
    bamfile : string
        The BAM file that needs processing

    ref_name : string
        The name of the reference sequence to collect the depths for

    min_depth : int
        The minimum depth to report coverage for (0 will be reported if coverage < min_depth at a given position)

    ignore_deletions : bool
        If true, positional depth counts will ignore reads with reference deletions

    Returns
    -------
    list
        Index is the reference position, value is the corresponding coverage depth
    """
    # check the BAM file exists
    if not os.path.exists(bamfile):
        raise Exception(f"bamfile doesn't exist {bamfile}")

    # open the BAM file
    bam_alignment = pysam.AlignmentFile(bamfile, "rb")

    # get the TID for the reference
    tid = bam_alignment.get_tid(ref_name)
    if tid == -1:
        raise Exception(f"bamfile does not contain specified reference {ref_name}")

    # create a depth vector to hold the depths at each reference position
    depths = [0] * bam_alignment.get_reference_length(ref_name)

    # generate the pileup
    for pileupcolumn in bam_alignment.pileup(
        ref_name,
        start=0,
        stop=bam_alignment.get_reference_length(ref_name),
        max_depth=100000000,
        truncate=False,
        min_base_quality=0,
    ):
        # process the pileup column
        for pileupread in pileupcolumn.pileups:

            # process the pileup read
            if pileupread.is_refskip:
                continue

            if pileupread.is_del:
                if not ignore_deletions:
                    depths[pileupcolumn.pos] += 1

            elif not pileupread.is_del:
                depths[pileupcolumn.pos] += 1

            else:
                raise Exception("unhandled pileup read encountered")

        # if final depth for pileup column < min_depth, report 0 and update the mask_vector
        if depths[pileupcolumn.pos] < min_depth:
            depths[pileupcolumn.pos] = 0

    # close file and return depth vector
    bam_alignment.close()
    return depths


# from https://www.geeksforgeeks.org/python-make-a-list-of-intervals-with-sequential-numbers/
def intervals_extract(iterable):
    iterable = sorted(set(iterable))
    for _, group in itertools.groupby(enumerate(iterable), lambda t: t[1] - t[0]):
        group = list(group)
        yield [group[0][1], group[-1][1]]


def go(args):

    # open the reference sequence and collect the sequence header and sequence length of the first record
    records = [x for x in SeqIO.parse(args.reference, "fasta")]
    intervals = []

    for record in records:
        seq_id = record.id
        seq_len = len(record.seq)

        # collect the depths from the pileup, replacing any depth<min_depth with 0
        depths = collect_depths(
            args.bamfile,
            seq_id,
            args.depth,
            args.ignore_deletions
        )

        # check the number of positions in the reported depths matches the reference sequence
        if len(depths) != seq_len:
            print("pileup length did not match expected reference sequence length")

        # create a mask_vector that records reference positions where depth < min_depth
        mask_vector = []
        for pos, depth in enumerate(depths):
            if depth == 0:
                mask_vector.append(pos)

        # get the intervals from the mask_vector
        intervals = list(intervals_extract(mask_vector))

        # create the mask outfile
        maskfh = open(args.outfile, "a")
        for i in intervals:
            maskfh.write(f"{seq_id}\t{i[0] + 1}\t{i[1] + 1}\n")
        maskfh.close()


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--depth", type=int, default=20)
    parser.add_argument(
        "--ignore-deletions",
        action="store_true",
        default=False,
        help="if set, positional depth counts will ignore reads with reference deletions (i.e. evaluates positional depths on ref matches, not read span",
    )
    parser.add_argument("reference")
    parser.add_argument("bamfile")
    parser.add_argument("outfile")

    args = parser.parse_args()
    go(args)


if __name__ == "__main__":
    main()
