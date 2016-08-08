#! /usr/bin/env python

import sys, subprocess

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

