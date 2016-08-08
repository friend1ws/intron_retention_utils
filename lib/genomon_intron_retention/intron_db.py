#! /usr/bin/env python

import gzip, subprocess

def generate_intron_retention_list(ref_gene_file, output_file, donor_size, acceptor_size, chr_name_list):

    ucsc2new_chr = {}
    if chr_name_list is not None:
        with open(chr_name_list, 'r') as hin:
            for line in hin:
                F = line.rstrip('\n').split('\t')
                ucsc2new_chr[F[0]] = F[1]

    donor_size_exon, donor_size_intron = [int(x) for x in donor_size.split(',')]
    acceptor_size_intron, acceptor_size_exon = [int(x) for x in acceptor_size.split(',')]

    key2junction, key2gene_id, key2exon_num = {}, {}, {}
    with gzip.open(ref_gene_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            chr = ucsc2new_chr[F[2]] if F[2] in ucsc2new_chr else F[2] 
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


