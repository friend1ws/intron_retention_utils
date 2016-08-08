#! /usr/bin/env python

import re, gzip
import pysam

def filterImproper(input_bam, output_bam, mapq_thres):

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
        if flags[1] != "1": continue

        # skip if either of the read pair is unmapped
        if flags[2] == "1" or flags[3] == "1": continue

        # skip supplementary alignment
        if flags[8] == "1" or flags[11] == "1": continue

        # skip duplicated reads
        if flags[10] == "1": continue

        # skip if below the minimum mapping quality
        if (read.mapq < mapq_thres): continue


        bamfile_out.write(read)



def summarize_edge(edge_bed, edge_broaden_bed, output_file, margin):

    edge2count = {}
    edge_broaden2count = {}

    with open(edge_bed, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[len(F) - 1] != "1": continue
            key = '\t'.join(F[12:18])
            if key not in edge2count: edge2count[key] = 0
            edge2count[key] = edge2count[key] + 1


    with open(edge_broaden_bed, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[len(F) - 1] != str(margin * 2): continue
            if F[17] == '+':
                F[13] = str(int(F[13]) + margin - 1)
                F[14] = str(int(F[14]) - margin)
            else:
                F[13] = str(int(F[13]) + margin)
                F[14] = str(int(F[14]) - margin + 1)

            key = '\t'.join(F[12:18])
            if key not in edge_broaden2count: edge_broaden2count[key] = 0
            edge_broaden2count[key] = edge_broaden2count[key] + 1


    hout = open(output_file, 'w')
    for key in list(set(edge2count.keys() + edge_broaden2count.keys())):
        count_edge = edge2count[key] if key in edge2count else 0
        count_edge_broaden = edge_broaden2count[key] if key in edge_broaden2count else 0
        print >> hout, key + '\t' + str(count_edge) + '\t' + str(count_edge_broaden)


def generate_intron_retention_list(ref_gene_file, output_file, donor_size, acceptor_size):

    donor_size_exon, donor_size_intron = [int(x) for x in donor_size.split(',')]
    acceptor_size_intron, acceptor_size_exon = [int(x) for x in acceptor_size.split(',')]

    key2junction, key2gene_id, key2exon_num = {}, {}, {}
    with gzip.open(ref_gene_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            chr = F[2]
            starts = [int(x) for x in F[9].split(',') if x != '']
            ends = [int(x) for x in F[10].split(',') if x != '']
            strand = F[3]
            exon_num = int(F[8])
            gene_id = F[1]
            symbol = F[12]

            
            for i in range(0, exon_num - 1):
                if strand == '+': # donor
                    key = '\t'.join([chr, str(ends[i] - donor_size_exon), str(ends[i] + donor_size_intron), symbol, "donor", strand])
                else: # acceptor
                    key = '\t'.join([chr, str(ends[i] - acceptor_size_exon), str(ends[i] + acceptor_size_intron), symbol, "acceptor", strand])
                
                junction = chr + ':' + str(ends[i]) + '-' + str(starts[i + 1] + 1)

                if key not in key2junction: key2junction[key] = []
                if key not in key2gene_id: key2gene_id[key] = []
                if key not in key2exon_num: key2exon_num[key] = []

                key2junction[key].append(junction)
                key2gene_id[key].append(gene_id)
                key2exon_num[key].append(str(i))


            for i in range(1, exon_num):
                if strand == '+': # acceptor 
                    key = '\t'.join([chr, str(starts[i] - acceptor_size_intron), str(starts[i] + acceptor_size_exon), symbol, "acceptor", strand])
                else: # donor 
                    key = '\t'.join([chr, str(starts[i] - donor_size_intron), str(starts[i] + donor_size_exon), symbol, "donor", strand])

                junction = chr + ':' + str(ends[i - 1]) + '-' + str(int(starts[i]) + 1)

                if key not in key2junction: key2junction[key] = []
                if key not in key2gene_id: key2gene_id[key] = []
                if key not in key2exon_num: key2exon_num[key] = []
                
                key2junction[key].append(junction) 
                key2gene_id[key].append(gene_id) 
                key2exon_num[key].append(str(i)) 


    hout = open(output_file + ".tmp", 'w')
    for key in sorted(key2junction):
        print >> hout, '\t'.join([key, ','.join(key2junction[key]), ','.join(key2gene_id[key]), ','.join(key2exon_num[key])])

    hout.close()

    hout = open(output_file, 'w')
    subprocess.call(["sort", "-k1,1", "-k2,2n", "-k3,3n", output_file + ".tmp"], stdout = hout)
    hout.close()

    subprocess.call(["rm", "-rf", output_file + ".tmp"])

