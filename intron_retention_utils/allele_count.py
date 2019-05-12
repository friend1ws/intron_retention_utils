#! /usr/bin/env python

from __future__ import print_function
import sys, re, pysam

from . import my_seq

mut_seq_margin = 100

# not UCSC bed format (-1 based start positions)
def generate_template_seq(output_file, reference, mut_chr, mut_start, mut_end, mut_ref, mut_alt, 
                          motif_chr, motif_start, motif_end, motif_type, motif_strand,
                          junc_list, donor_size, acceptor_size, template_size):

    donor_size_exon, donor_size_intron = [int(x) for x in donor_size.split(',')]
    acceptor_size_intron, acceptor_size_exon = [int(x) for x in acceptor_size.split(',')]

    unique_junc_list = list(set(junc_list.split(',')))

    if motif_type == "donor":
        if motif_strand == '+':
            motif_exon_start, motif_exon_end = motif_start, motif_start + donor_size_exon - 1
            motif_intron_start, motif_intron_end = motif_start + donor_size_exon, motif_start + donor_size_exon + donor_size_intron - 1
        else:
            motif_intron_start, motif_intron_end = motif_start, motif_start + donor_size_intron - 1
            motif_exon_start, motif_exon_end = motif_start + donor_size_intron, motif_start + donor_size_exon + donor_size_intron - 1
    else: # acceptor
        if motif_strand == '+':
            motif_intron_start, motif_intron_end = motif_start, motif_start + acceptor_size_intron - 1
            motif_exon_start, motif_exon_end = motif_start + acceptor_size_intron, motif_start + acceptor_size_exon + acceptor_size_intron - 1
        else:
            motif_exon_start, motif_exon_end = motif_start, motif_start + acceptor_size_exon - 1
            motif_intron_start, motif_intron_end = motif_start + acceptor_size_exon, motif_start + acceptor_size_exon + acceptor_size_intron - 1

    key2seq = {}  
    # annotated splice junction without mutation
    cnum = 0
    for junc in sorted(unique_junc_list):
        junc_match = re.match(r'([^ \t\n\r\f\v,]+)\:(\d+)\-(\d+)', junc)
        junc_chr, junc_start, junc_end = junc_match.group(1), int(junc_match.group(2)), int(junc_match.group(3))

        seq = my_seq.get_seq(reference, junc_chr, junc_start - template_size + 1, junc_start) + \
              my_seq.get_seq(reference, junc_chr, junc_end, junc_end + template_size - 1)
        key2seq["splice_junction_negative_" + str(cnum)] = seq
        cnum = cnum + 1


    # annotated splice junction with mutation (only when mutations occur within exonic motif region)
    if mut_start >= motif_exon_start and mut_end <= motif_exon_end:
        cnum = 0
        for junc in sorted(unique_junc_list):
            junc_match = re.match(r'([^ \t\n\r\f\v,]+)\:(\d+)\-(\d+)', junc)
            junc_chr, junc_start, junc_end = junc_match.group(1), int(junc_match.group(2)), int(junc_match.group(3))

            if (motif_type == "donor" and motif_strand == '+') or (motif_type == "acceptor" and motif_strand == '-'):
                mut_seq_tmp = my_seq.get_seq(reference, junc_chr, junc_start - mut_seq_margin + 1, junc_start)
                mut_start_rel, mut_end_rel = mut_start - junc_start + mut_seq_margin - 1, mut_end - junc_start + mut_seq_margin - 1
            else:
                mut_seq_tmp = my_seq.get_seq(reference, junc_chr, junc_end, junc_end + mut_seq_margin - 1)
                mut_start_rel, mut_end_rel = mut_start - junc_end, mut_end - junc_end

            # for debug
            if mut_ref != '-' and mut_seq_tmp[mut_start_rel:(mut_end_rel + 1)] != mut_ref:
                print('\t'.join([mut_chr, str(mut_start), str(mut_end), mut_ref, mut_alt, junc]), file = sys.stderr)
                print('\t'.join([mut_seq_tmp[mut_start_rel:(mut_end_rel + 1)], mut_ref]), file = sys.stderr)
                print("mutation inconsistent!!!", file = sys.stderr)
                sys.exit(1) 

            # SNV
            if mut_ref != '-' and mut_alt != '-': mut_seq_tmp = mut_seq_tmp[:mut_start_rel] + mut_alt + mut_seq_tmp[(mut_end_rel + 1):]

            # deletion
            if mut_alt == '-': mut_seq_tmp = mut_seq_tmp[:mut_start_rel] + mut_seq_tmp[(mut_end_rel + 1):]

            # insertion
            if mut_ref == '-': mut_seq_tmp = mut_seq_tmp[:(mut_start_rel + 1)] + mut_alt + mut_seq_tmp[(mut_start_rel + 1):]

            if (motif_type == "donor" and motif_strand == '+') or (motif_type == "acceptor" and motif_strand == '-'):
                seq = mut_seq_tmp[(-template_size):] + my_seq.get_seq(reference, junc_chr, junc_end, junc_end + template_size - 1)
            else:
                seq = my_seq.get_seq(reference, junc_chr, junc_start - template_size + 1, junc_start) + mut_seq_tmp[:(template_size)]

            key2seq["splice_junction_positive_" + str(cnum)] = seq
            cnum = cnum + 1

            
    # intron retention without mutation
    if (motif_type == "donor" and motif_strand == '+') or (motif_type == "acceptor" and motif_strand == '-'):
        seq_left_tmp = my_seq.get_seq(reference, junc_chr, motif_exon_end - mut_seq_margin + 1, motif_exon_end)
        seq_right_tmp = my_seq.get_seq(reference, junc_chr, motif_intron_start, motif_intron_start + mut_seq_margin - 1)
        boundary_pos = motif_exon_end
    else:
        seq_left_tmp = my_seq.get_seq(reference, junc_chr, motif_intron_end - mut_seq_margin + 1, motif_intron_end)
        seq_right_tmp = my_seq.get_seq(reference, junc_chr, motif_exon_start, motif_exon_start + mut_seq_margin - 1)
        boundary_pos = motif_intron_end
    key2seq["intron_retention_negative"] = seq_left_tmp[(-template_size):] + seq_right_tmp[:(template_size)]

   
    # intron retention with mutation
    mut_seq_left_tmp, mut_seq_right_tmp = seq_left_tmp, seq_right_tmp

    # in this case, we remove nucleotides from concatenated sequences
    mut_start_rel, mut_end_rel = mut_start - boundary_pos + mut_seq_margin - 1, mut_end - boundary_pos + mut_seq_margin - 1
    mut_seq_tmp = mut_seq_left_tmp + mut_seq_right_tmp


    # SNV
    if mut_ref != '-' and mut_alt != '-': 

        # for debug
        if mut_seq_tmp[mut_start_rel] != mut_ref:
            print('\t'.join([mut_chr, str(mut_start), str(mut_end), mut_ref, mut_alt]))
            print(mut_seq_tmp)
            print(mut_start_rel)
            print("mutation inconsistent!!!", file = sys.stderr)
            sys.exit(1)
        mut_seq_tmp = mut_seq_tmp[:mut_start_rel] + mut_alt + mut_seq_tmp[(mut_end_rel + 1):]

        mut_seq_start_pos = mut_seq_margin - template_size
        mut_seq_end_pos = mut_seq_start_pos + 2 * template_size
        key2seq["intron_retention_positive"] = mut_seq_tmp[mut_seq_start_pos:mut_seq_end_pos]


    elif mut_alt == '-': # deletion

        # for debug
        if mut_seq_tmp[mut_start_rel:(mut_end_rel + 1)] != mut_ref != '-':
            print('\t'.join([mut_chr, str(mut_start), str(mut_end), mut_ref, mut_alt]))
            print("mutation inconsistent!!!", file = sys.stderr)
            sys.exit(1)
        mut_seq_tmp = mut_seq_tmp[:mut_start_rel] + mut_seq_tmp[(mut_end_rel + 1):]
 
        del_size_left = max(0, (min(boundary_pos, mut_end) - mut_start))
        mut_seq_start_pos = mut_seq_margin - template_size - del_size_left
        mut_seq_end_pos = mut_seq_start_pos + 2 * template_size

        key2seq["intron_retention_positive"] = mut_seq_tmp[mut_seq_start_pos:mut_seq_end_pos]


    elif mut_ref == '-': #insertion
     
        mut_seq_tmp = mut_seq_tmp[:(mut_start_rel + 1)] + mut_alt + mut_seq_tmp[(mut_start_rel + 1):]
        ins_size_left = len(mut_alt)
        mut_seq_start_pos = mut_seq_margin - template_size + ins_size_left
        mut_seq_end_pos = mut_seq_start_pos + 2 * template_size
        key2seq["intron_retention_positive"] = mut_seq_tmp[mut_seq_start_pos:mut_seq_end_pos]


    hout = open(output_file, 'w')
    for key in sorted(key2seq):
        print(">" + key + '\n' + key2seq[key], file = hout)
    
    hout.close()
 

    
def extract_read_around_boundary(bam_file, output_file, reference, motif_chr, motif_start, motif_end, read_search_margin):

    bamfile = pysam.Samfile(bam_file, 'rb')
    hout = open(output_file, 'w') 
    for read in bamfile.fetch(motif_chr, max(0, motif_start - read_search_margin), motif_end + read_search_margin):

        # get the flag information
        flags = format(int(read.flag), "#014b")[:1:-1]

        # skip unmapped read 
        if flags[2] == "1" or flags[3] == "1": continue 

        # skip supplementary alignment
        if flags[8] == "1" or flags[11] == "1": continue

        # skip duplicated reads
        if flags[10] == "1": continue

        read_id = read.qname + '/1' if flags[6] == '1' else read.qname + '/2'
        
        print('>' + read_id + '\n' + read.seq, file = hout)

    bamfile.close()
    hout.close()
 
