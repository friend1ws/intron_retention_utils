#! /usr/bin/env python

from __future__ import print_function
import re, gzip
import pysam

def filterImproper(input_bam, output_bam, keep_improper_pair):

    """
    This function is used for filtering short read with improper pairs,
    and low mapping qualities.
    """

    bamfile_in = pysam.AlignmentFile(input_bam, 'rb')
    bamfile_out = pysam.AlignmentFile(output_bam, 'wb', template = bamfile_in)

    for read in bamfile_in.fetch():
    
        # get the flag information
        flags = format(int(read.flag), "#014b")[:1:-1]

        # skip improper read pair
        if flags[1] != "1" and not keep_improper_pair: continue

        # skip if either of the read pair is unmapped
        if flags[2] == "1" or flags[3] == "1": continue

        # skip supplementary alignment
        if flags[8] == "1" or flags[11] == "1": continue

        # skip duplicated reads
        if flags[10] == "1": continue

        bamfile_out.write(read)



def summarize_edge(edge_bed, edge_broaden_bed, output_file, margin, mapq_thres):

    edge2count = {}
    edge_broaden2count = {}

    """
    with open(edge_bed, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if int(F[4]) < mapq_thres: continue 
            if F[len(F) - 1] != "1": continue
            key = F[12] + '\t' + '\t'.join(F[14:(len(F) - 1)])
            if key not in edge2count: edge2count[key] = 0
            edge2count[key] = edge2count[key] + 1
    """

    with open(edge_broaden_bed, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if int(F[4]) < mapq_thres: continue
            # if F[len(F) - 1] != str(margin * 2): continue
            if (F[16] == "donor" and F[17] == '+') or (F[16] == "acceptor" and F[17] == '-'): 
                F[13] = str(int(F[13]) + margin - 1)
                F[14] = str(int(F[14]) - margin)
            else:
                F[13] = str(int(F[13]) + margin)
                F[14] = str(int(F[14]) - margin + 1)

            key = F[12] + '\t' + '\t'.join(F[14:(len(F) - 1)])
            if F[len(F) - 1] == str(margin * 2):
                if key not in edge_broaden2count: edge_broaden2count[key] = 0
                edge_broaden2count[key] = edge_broaden2count[key] + 1

    
            # print(F[10], F[11])
            block_sizes = [int(x) for x in F[10].split(',') if x != '']
            block_starts = [int(x) for x in F[11].split(',') if x != '']

            for i in range(len(block_sizes)):
                # if block_starts[i] == '': continue
                tstart = int(F[6]) + block_starts[i] + 1
                tend = int(F[6]) + block_starts[i] + block_sizes[i]
                if tstart <= int(F[14]) and int(F[14]) <= tend:
                    if key not in edge2count: edge2count[key] = 0
                    edge2count[key] = edge2count[key] + 1

    hout = open(output_file, 'w')
    for key in list(set(list(edge2count) + list(edge_broaden2count))):
        count_edge = edge2count[key] if key in edge2count else 0
        count_edge_broaden = edge_broaden2count[key] if key in edge_broaden2count else 0
        print(key + '\t' + str(count_edge) + '\t' + str(count_edge_broaden), file = hout)

    hout.close()

