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


def process_result(input_file, intron_retention_header_file, output_file):

    hout = open(output_file, 'w')

    with open(intron_retention_header_file, 'r') as hin:
        intron_retention_header = hin.readline().rstrip('\n')
        

    print >> hout, '\t'.join(["Chr_Mut", "Start_Mut", "End_Mut", "Ref_Mut", "Alt_Mut", "Intron_Retention_Type"]) + '\t' + intron_retention_header

    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            print >> hout, '\t'.join(F[0:5]) + '\t' + '\t'.join(F[8:])

    hout.close()

    
