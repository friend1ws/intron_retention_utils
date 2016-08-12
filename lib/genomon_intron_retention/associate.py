#! /usr/bin/env python

import sys, re, subprocess

def generate_mutation_target(input_file, output_file, output_header, donor_size, acceptor_size):

    donor_size_exon, donor_size_intron = [int(x) for x in donor_size.split(',')]
    acceptor_size_intron, acceptor_size_exon = [int(x) for x in acceptor_size.split(',')]
    
    header2ind = {}
    hout = open(output_file + ".tmp", 'w')
    with open(input_file, 'r') as hin:
        
        header = hin.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            header2ind[cname] = i

        hout_header = open(output_header, 'w')
        print >> hout_header, '\t'.join(header)
        hout_header.close()


        for line in hin:
            F = line.rstrip('\n').split('\t')

            boundary_pos = int(F[header2ind["Boundary_Pos"]])
            motif_type = F[header2ind["Motif_Type"]]
            strand = F[header2ind["Strand"]]

            # direct position
            if motif_type == "donor":
                if strand == '+':
                    start, end = boundary_pos - donor_size_exon, boundary_pos + donor_size_intron
                else:
                    start, end = boundary_pos - donor_size_intron - 1, boundary_pos + donor_size_exon - 1
            else:
                if strand == '+':
                    start, end = boundary_pos - acceptor_size_intron - 1, boundary_pos + acceptor_size_exon - 1
                else:
                    start, end = boundary_pos - acceptor_size_exon, boundary_pos + acceptor_size_intron

            print >> hout, '\t'.join([F[header2ind["Chr"]], str(start), str(end), "direct-impact"] + F)

            # opposite-side-impact 
            unique_junc_list = list(set(F[header2ind["Junction_List"]].split(',')))
            for junc in sorted(unique_junc_list):

                junc_match = re.match(r'([^ \t\n\r\f\v,]+)\:(\d+)\-(\d+)', junc)
                if junc_match is None: 
                    print >> sys.stderr, "junction mismatch at " + '\n' + '\t'.join(F)
                    sys.exit(1)

                opposite_motif_type = "donor" if motif_type == "acceptor" else "acceptor"

                if (motif_type == "donor" and strand == '+') or (motif_type == "acceptor" and strand == '-'):
                    opposite_boundary_pos = int(junc_match.group(3))
                else:
                    opposite_boundary_pos = int(junc_match.group(2)) 

                if opposite_motif_type == "donor":
                    if strand == '+':
                        start, end = opposite_boundary_pos - donor_size_exon, opposite_boundary_pos + donor_size_intron
                    else:
                        start, end = opposite_boundary_pos - donor_size_intron - 1, opposite_boundary_pos + donor_size_exon - 1
                else:   
                    if strand == '+':
                        start, end = opposite_boundary_pos - acceptor_size_intron - 1, opposite_boundary_pos + acceptor_size_exon - 1
                    else:
                        start, end = opposite_boundary_pos - acceptor_size_exon, opposite_boundary_pos + acceptor_size_intron
                        
                print >> hout, '\t'.join([F[header2ind["Chr"]], str(start), str(end), "opposite-side-impact"] + F)
 
    hout.close()

    hout = open(output_file, 'w')
    subprocess.call(["sort", "-k1,1", "-k2,2n", "-k3,3n", output_file + ".tmp"], stdout = hout)
    hout.close()

    subprocess.call(["rm", "-rf", output_file + ".tmp"])


def process_result(input_file, intron_retention_header_file, output_file, donor_size, acceptor_size):

    donor_size_exon, donor_size_intron = [int(x) for x in donor_size.split(',')]
    acceptor_size_intron, acceptor_size_exon = [int(x) for x in acceptor_size.split(',')]


    hout = open(output_file, 'w')

    with open(intron_retention_header_file, 'r') as hin:
        intron_retention_header = hin.readline().rstrip('\n')
        

    print >> hout, intron_retention_header + '\t' + '\t'.join(["Mutation_Key", "Motif_Pos", "Mutation_Type", "Is_Canonical", "Intron_Retention_Type"])

    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            mutation_key = F[3]
            motif_pos = F[4] + ':' + str(int(F[5]) + 1) + '-' + F[6] + ',' + F[12]
            
            if (F[11] == "donor" and F[7] == "direct-impact") or (F[11] == "acceptor" and F[7] == "opposite-side-impact"):
                motif_type = "splicing donor disruption"
            elif (F[11] == "acceptor" and F[7] == "direct-impact") or (F[11] == "donor" and F[7] == "opposite-side-impact"): 
                motif_type = "splicing acceptor disruption"

            # check cannonical
            mut_start, mut_end = int(F[1]) + 1, int(F[2])
            if motif_type == "splicing donor disruption":
                if F[12] == '+':
                    canonical_start, canonical_end = int(F[5]) + donor_size_exon + 1, int(F[5]) + donor_size_exon + 2
                else:
                    canonical_start, canonical_end = int(F[5]) + donor_size_intron - 1, int(F[5]) + donor_size_intron 
            else:
                if F[12] == '+':
                    canonical_start, canonical_end = int(F[5]) + acceptor_size_intron - 1, int(F[5]) + acceptor_size_intron
                else:
                    canonical_start, canonical_end = int(F[5]) + acceptor_size_exon + 1, int(F[5]) + acceptor_size_exon + 2

            is_canonical = "canonical" if (mut_start <= canonical_end and mut_end >= canonical_start) else "non-canonical" 

            print >> hout, '\t'.join(F[8:]) + '\t' + mutation_key + '\t' + motif_pos + '\t' + motif_type + '\t' + is_canonical + '\t' + F[7]

    hout.close()

    
