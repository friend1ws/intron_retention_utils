#! /usr/bin/env python

import sys, gzip

junction2annot = {}


def proc_line(line, type):

    F  = line.rstrip('\n').split('\t')
    
    chr = F[2]
    starts = F[9].split(',')
    ends = F[10].split(',')
    strand = F[3]
    exon_num = int(F[8])
    gene = F[1]
    symbol = F[12]
    
    for i in range(1, exon_num):
        
        key = chr + '\t' + starts[i] + '\t' + str(int(starts[i]) + 1) + '\t' + '-'
        annot = symbol + '(' + gene + ').' if type == "refGene" else gene
        annot = annot + (str(i) + ".start" if strand == '+' else str(exon_num - i) + ".end")
        if key in junction2annot: 
            junction2annot[key] = junction2annot[key] + ',' + annot
        else:
            junction2annot[key] = annot
    
    
    for i in range(0, exon_num - 1):
        
        key = chr + '\t' + str(int(ends[i]) - 1) + '\t' + ends[i] + '\t' + '+' 
        annot = symbol + '(' + gene + ').' if type == "refGene" else gene
        annot = symbol + '(' + gene + ').' + (str(i) + ".end" if strand == '+' else str(exon_num - i - 1) + ".start")
        
        if key in junction2annot: 
            junction2annot[key] = junction2annot[key] + ',' + annot
        else:
            junction2annot[key] = annot


print >> sys.stderr, "Start reading refGene.txt.gz"
num = 0
with gzip.open("refGene.txt.gz", 'r') as hin:
    for line in hin:
        proc_line(line, "refGene")
        num = num + 1
        if num % 1000 == 0: print >> sys.stderr, str(num) + " genes completed."


print >> sys.stderr, "Start reading ensGene.txt.gz"
num = 0
with gzip.open("ensGene.txt.gz", 'r') as hin:
    for line in hin:
        proc_line(line, "ensGene")
        num = num + 1
        if num % 1000 == 0: print >> sys.stderr, str(num) + " genes completed."


for junction in sorted(junction2annot):
    junc_info = junction.split('\t')
    print '\t'.join(junc_info[0:3]) + '\t' + junction2annot[junction] + '\t' + '0' + '\t' + junc_info[3]
