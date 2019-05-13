#! /usr/bin/env python

from __future__ import print_function
import sys, subprocess
import pysam

from . import my_seq

def anno2vcf(input_file, output_file, reference):

    hin = open(input_file, 'r')
    hout = open(output_file, 'w')

    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[0].startswith('#'): continue
        if F[0] == "Chr" and F[1] == "Start" and F[2] == "End" and F[3] == "Ref" and F[4] == "Alt": continue

        pos, ref, alt = F[1], F[3], F[4]
        
        # insertion
        if F[3] == "-":
            # get the sequence for the reference base
            seq = my_seq.get_seq(reference, F[0], int(F[1]), int(F[1]))
            ref, alt = seq, seq + F[4]

        # deletion
        if F[4] == "-":
            # get the sequence for the reference base
            seq = my_seq.get_seq(reference, F[0], int(F[1]) - 1, int(F[1]) - 1)
            pos, ref, alt = str(int(F[1]) - 1), seq + F[3], seq

        # QUAL = int(float(F[15]) * 10)
        QUAL = 60
        # INFO = "TD=" + F[9] + ";TV=" + F[10] + ";ND=" + F[13] + ";NV=" + F[14] + ";SOMATIC"
        INFO = "SOMATIC"

        print(F[0] + "\t" + pos + "\t.\t" + ref + "\t" + alt \
           + "\t" + str(QUAL) + "\t" + "PASS" + "\t" + INFO, file = hout)
    
    hin.close()
    hout.close()


def remove_vcf_header(input_file, output_file):

    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')
            if not line.startswith('#'): print(line, file = hout)

    hout.close()


def vcf2bed(mutation_file, output_file):

    hout = open(output_file, 'w')
    with open(mutation_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            tchr, tstart, tend = F[0], str(int(F[1]) - 1), F[1]

            # deletion
            if len(F[3]) > 1:
                tstart, tend = F[1], str(int(F[1]) + len(F[3]) - 1)
            # insertion
            elif len(F[4]) > 1:
                tstart, tend = str(int(F[1]) - 1), F[1]
            # substitution
            elif len(F[3]) == 1 and len(F[4]) == 1:
                tstart, tend = str(int(F[1]) - 1), F[1]
 
            print(tchr + '\t' + tstart + '\t' + tend + '\t' + ','.join([F[0], F[1], F[3], F[4]]), file = hout)

    hout.close()


def anno2bed(mutation_file, output_file):

    hout = open(output_file + ".tmp", 'w')
    with open(mutation_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            if F[0].startswith('#'): continue
            if F[0] == "Chr" and F[1] == "Start" and F[2] == "End" and F[3] == "Ref" and F[4] == "Alt": continue

            print('\t'.join([F[0], str(int(F[1]) - 1), F[2], F[3], F[4]]), file = hout)

    hout.close()

    hout = open(output_file, 'w')
    subprocess.check_call(["sort", "-k1,1", "-k2,2n", "-k3,3n", output_file + ".tmp"], stdout = hout)
    hout.close()

    subprocess.check_call(["rm", "-rf", output_file + ".tmp"])


def genosv2bed(input_file, output_file):

    hout = open(output_file, 'w')
    num = 1
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0].startswith('#'): continue
            if F[0] == "Chr_1" and F[1] == "Pos_1": continue
            chr1, chr2 = F[0], F[3]
            start1, end1 = str(int(F[1]) - 1), F[1]
            start2, end2 = str(int(F[4]) - 1), F[4]
            dir1, dir2 = F[2], F[5]
            inseq = F[6]

            key = ','.join([chr1, end1, dir1, chr2, end2, dir2, inseq])

            print('\t'.join([chr1, start1, end1, dir1, key]), file = hout)
            print('\t'.join([chr2, start2, end2, dir2, key]), file = hout)

    hout.close()


