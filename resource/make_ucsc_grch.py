#! /usr/bin/env python

import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

hout = open(output_file, 'w') 
with open(input_file, 'r') as hin:
    for line in hin:
        if line.startswith('#'): continue
        F = line.rstrip('\n\r').split('\t')
        if F[4].startswith('CM'):
            print >> hout, F[9] + '\t' + F[2]
        else:
            print >> hout, F[9] + '\t' + F[4]

hout.close()
