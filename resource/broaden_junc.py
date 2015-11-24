#! /usr/bin/env python

import sys

input_file = sys.argv[1]
margin = int(sys.argv[2])

with open(input_file, 'r') as hin:
    for line in hin:
        F  = line.rstrip('\n').split('\t')

        if F[5] == '+':
            print F[0] + '\t' + str(int(F[1]) - margin + 1) + '\t' + str(int(F[2]) + margin) + '\t' + '\t'.join(F[3:])
        else:
            print F[0] + '\t' + str(int(F[1]) - margin) + '\t' + str(int(F[2]) + margin - 1) + '\t' + '\t'.join(F[3:])


