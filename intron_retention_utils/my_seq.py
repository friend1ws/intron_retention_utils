#! /usr/bin/env python

from __future__ import print_function
import sys, re
import pysam

def get_seq(reference, chr, start, end):

    seq = ""
    for item in pysam.faidx(reference, chr + ":" + str(start) + "-" + str(end)):
        # if item[0] == ">": continue
        seq = seq + item.rstrip('\n')
    seq = seq.replace('>', '')
    seq = seq.replace(chr + ":" + str(start) + "-" + str(end), '')

    if re.search(r'[^ACGTNacgtn]', seq) is not None:
        print("The return value in get_seq function includes non-nucleotide characters:", file = sys.stderr)
        print(seq, file = sys.stderr)
        sys.exit(1)


    return seq


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
                  'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}

    return("".join(complement.get(base, base) for base in reversed(seq)))


