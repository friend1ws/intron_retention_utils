#! /usr/bin/env python

import sys, subprocess
import pysam

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
            seq = ""    
            for item in pysam.faidx(reference, F[0] + ":" + str(F[1]) + "-" + str(F[1])):
                if item[0] == ">": continue
                seq = seq + item.rstrip('\n')
            ref, alt = seq, seq + F[4]

        # deletion
        if F[4] == "-":
            # get the sequence for the reference base
            seq = ""    
            for item in pysam.faidx(reference, F[0] + ":" + str(int(F[1]) - 1) + "-" + str(int(F[1]) - 1)):
                if item[0] == ">": continue
                seq = seq + item.rstrip('\n')

            pos, ref, alt = str(int(F[1]) - 1), seq + F[3], seq

        # QUAL = int(float(F[15]) * 10)
        QUAL = 60
        # INFO = "TD=" + F[9] + ";TV=" + F[10] + ";ND=" + F[13] + ";NV=" + F[14] + ";SOMATIC"
        INFO = "SOMATIC"

        print >> hout, F[0] + "\t" + pos + "\t.\t" + ref + "\t" + alt \
            + "\t" + str(QUAL) + "\t" + "PASS" + "\t" + INFO 


def remove_vcf_header(input_file, output_file):

    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')
            if not line.startswith('#'): print >> hout, line

    hout.close()


def vcf2bed(mutation_file, output_file):

    hout = open(output_file, 'w')
    with open(mutation_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            tchr, tstart, tend = F[0], str(int(F[1]) - 1), F[1]

            # deletion
            if len(F[3]) > 1:
                tstart, tend = F[1], str(int(F[1]) + len(F[3]))
            # insertion
            elif len(F[4]) > 1:
                tstart, tend = str(int(F[1]) - 1), F[1]
            # substitution
            elif len(F[3]) == 1 and len(F[4]) == 1:
                tstart, tend = str(int(F[1]) - 1), F[1]
 
            print >> hout, tchr + '\t' + tstart + '\t' + tend + '\t' + ','.join([F[0], F[1], F[3], F[4]])

    hout.close()


def anno2bed(mutation_file, output_file):

    hout = open(output_file + ".tmp", 'w')
    with open(mutation_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            if F[0].startswith('#'): continue
            if F[0] == "Chr" and F[1] == "Start" and F[2] == "End" and F[3] == "Ref" and F[4] == "Alt": continue

            print >> hout, '\t'.join([F[0], str(int(F[1]) - 1), F[2], F[3], F[4]])

    hout.close()

    hout = open(output_file, 'w')
    subprocess.call(["sort", "-k1,1", "-k2,2n", "-k3,3n", output_file + ".tmp"], stdout = hout)
    hout.close()

    subprocess.call(["rm", "-rf", output_file + ".tmp"])

