#! /usr/bin/env python

from __future__ import print_function
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
        print('\t'.join(header), file = hout_header)
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

            print('\t'.join([F[header2ind["Chr"]], str(start), str(end), "Direct impact"] + F), file = hout)

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
                        
                print('\t'.join([F[header2ind["Chr"]], str(start), str(end), "Opposite side impact"] + F), file = hout)
 
    hout.close()

    hout = open(output_file, 'w')
    subprocess.check_call(["sort", "-f", "-k1,1", "-k2,2n", "-k3,3n", output_file + ".tmp"], stdout = hout)
    hout.close()

    subprocess.check_call(["rm", "-rf", output_file + ".tmp"])


def generate_sv_target(input_file, output_file, output_header, intron_margin):

    header2ind = {}
    hout = open(output_file + ".tmp", 'w')
    with open(input_file, 'r') as hin:

        header = hin.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            header2ind[cname] = i

        hout_header = open(output_header, 'w')
        print('\t'.join(header), file = hout_header)
        hout_header.close()

        for line in hin:
            F = line.rstrip('\n').split('\t')

            unique_junc_list = list(set(F[header2ind["Junction_List"]].split(',')))
            for junc in sorted(unique_junc_list):

                junc_match = re.match(r'([^ \t\n\r\f\v,]+)\:(\d+)\-(\d+)', junc)
                if junc_match is None:
                    print("junction mismatch at " + '\n' + '\t'.join(F), file = sys.stderr)
                    sys.exit(1)

                jchr, jstart, jend = junc_match.group(1), int(junc_match.group(2)), int(junc_match.group(3))

                print('\t'.join([jchr, str(jstart - intron_margin), str(jend + intron_margin)] + F), file = hout)

    hout.close()


    hout = open(output_file, 'w')
    subprocess.check_call(["sort", "-f", "-k1,1", "-k2,2n", "-k3,3n", output_file + ".tmp"], stdout = hout)
    hout.close()

    subprocess.check_call(["rm", "-rf", output_file + ".tmp"])


def process_result(input_file, intron_retention_header_file, output_file, donor_size, acceptor_size):

    donor_size_exon, donor_size_intron = [int(x) for x in donor_size.split(',')]
    acceptor_size_intron, acceptor_size_exon = [int(x) for x in acceptor_size.split(',')]


    hout = open(output_file, 'w')

    with open(intron_retention_header_file, 'r') as hin:
        intron_retention_header = hin.readline().rstrip('\n')
        

    print(intron_retention_header + '\t' + '\t'.join(["Mutation_Key", "Motif_Pos", "Mutation_Type", "Is_Canonical", "Intron_Retention_Type"]), file = hout)

    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            mutation_key = F[3]
            motif_pos = F[4] + ':' + str(int(F[5]) + 1) + '-' + F[6] + ',' + F[12]
            
            if (F[11] == "donor" and F[7] == "Direct impact") or (F[11] == "acceptor" and F[7] == "Opposite side impact"):
                motif_type = "Donor disruption"
            elif (F[11] == "acceptor" and F[7] == "Direct impact") or (F[11] == "donor" and F[7] == "Opposite side impact"): 
                motif_type = "Acceptor disruption"

            # check cannonical
            mut_start, mut_end = int(F[1]) + 1, int(F[2])
            if motif_type == "Donor disruption":
                if F[12] == '+':
                    canonical_start, canonical_end = int(F[5]) + donor_size_exon + 1, int(F[5]) + donor_size_exon + 2
                else:
                    canonical_start, canonical_end = int(F[5]) + donor_size_intron - 1, int(F[5]) + donor_size_intron 
            else:
                if F[12] == '+':
                    canonical_start, canonical_end = int(F[5]) + acceptor_size_intron - 1, int(F[5]) + acceptor_size_intron
                else:
                    canonical_start, canonical_end = int(F[5]) + acceptor_size_exon + 1, int(F[5]) + acceptor_size_exon + 2

            is_canonical = "Canonical" if (mut_start <= canonical_end and mut_end >= canonical_start) else "Non-canonical" 

            print('\t'.join(F[8:]) + '\t' + mutation_key + '\t' + motif_pos + '\t' + motif_type + '\t' + is_canonical + '\t' + F[7], file = hout)

    hout.close()



def process_result_sv(input_file, intron_retention_header_file, output_file):

    hout = open(output_file, 'w')

    with open(intron_retention_header_file, 'r') as hin:
        intron_retention_header = hin.readline().rstrip('\n')


    print(intron_retention_header + '\t' + "SV_Key", file = hout)

    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            sv_key = F[4]

            """
            # check direction consistency
            consistent_flag = 0

            if F[11] not in ["donor", "acceptor"]: print "error!!"
            if F[12] not in ['+', '-']: print "error!!"

            if F[11] == "donor" and F[12] == "+" and F[3] == "-": consistent_flag = 1
            if F[11] == "donor" and F[12] == "-" and F[3] == "+": consistent_flag = 1
            if F[11] == "acceptor" and F[12] == "+" and F[3] == "+": consistent_flag = 1
            if F[11] == "acceptor" and F[12] == "-" and F[3] == "-": consistent_flag = 1

            if consistent_flag == 1: continue
            """
            print('\t'.join(F[8:]) + '\t' + sv_key, file = hout)

    hout.close()

