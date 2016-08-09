#! /usr/bin/env python

import sys, os, subprocess
import utils, intron_db, mutation, allele_count, pyssw

def simple_count_main(args):

    input_bam = args.bam_file
    output_prefix = args.output_prefix
    annotation_dir = args.annotation_dir
    bedtools_path = args.bedtools_path
    mapq_thres = args.q

    output_prefix_dir = os.path.dirname(output_prefix)
    if output_prefix_dir != "" and not os.path.exists(output_prefix_dir):
       os.makedirs(output_prefix_dir)


    utils.filterImproper(input_bam, output_prefix + ".filt.bam", mapq_thres)


    hOUT = open(output_prefix + ".filt.bed12", 'w')
    s_ret = subprocess.call([bedtools_path + "/bedtools", "bamtobed", "-bed12", "-i", output_prefix + ".filt.bam"], stdout = hOUT)
    hOUT.close()

    if s_ret != 0:
        print >> sys.stderr, "error in generating filt.bed12"
        sys.exit(1)


    hOUT = open(output_prefix + ".edge.bed", 'w')
    s_ret = subprocess.call([bedtools_path + "/bedtools", "intersect", "-a", output_prefix + ".filt.bed12", "-b", annotation_dir + "/edge.bed", 
                     "-split", "-wao"], stdout = hOUT)
    hOUT.close()

    if s_ret != 0:
        print >> sys.stderr, "error in generating edge.bed"
        sys.exit(1)

    hOUT = open(output_prefix + ".edge_broaden.bed", 'w')
    s_ret = subprocess.call([bedtools_path + "/bedtools", "intersect", "-a", output_prefix + ".filt.bed12", "-b", annotation_dir + "/edge_broaden.bed", 
                     "-split", "-wao"], stdout = hOUT)
    hOUT.close()

    if s_ret != 0:
        print >> sys.stderr, "error in generating edge_broaden.bed"
        sys.exit(1)

    utils.summarize_edge(output_prefix + ".edge.bed", output_prefix + ".edge_broaden.bed", output_prefix + ".intron.bed", 5)


    subprocess.call(["rm", "-rf", output_prefix + ".filt.bam"])
    subprocess.call(["rm", "-rf", output_prefix + ".filt.bed12"])
    subprocess.call(["rm", "-rf", output_prefix + ".exon.bed"])
    subprocess.call(["rm", "-rf", output_prefix + ".exon2base.txt"])



def allele_count_main(args):

    output_dir = os.path.dirname(args.output_file)
    if output_dir != "" and not os.path.exists(output_dir):
       os.makedirs(output_dir)

    intron_db.generate_intron_retention_list(args.ref_gene_file, args.output_file + ".intron_retention_list.bed", 
                                             args.donor_size, args.acceptor_size, args.chr_name_list)

    mutation.anno2bed(args.mutation_file, args.output_file + ".mutation_list.bed")

    hout = open(args.output_file + ".mutation_list.overlap.bed", 'w')
    subprocess.call(["bedtools", "intersect", "-a", args.output_file + ".mutation_list.bed",
                     "-b", args.output_file + ".intron_retention_list.bed", "-wa", "-wb"], stdout = hout)
    hout.close()
    
    cnum = 0
    hout = open(args.output_file, 'w')
    # print header
    print >> hout, '\t'.join(["Gene_Symbol", "Chr_Mut", "Start_Mut", "End_Mut", "Ref_Mut", "Alt_Mut", 
                              "Chr_Motif", "Start_Motif", "End_Motif", "Type_Motif", "Strand_Motif",
                              "Splice_Junction_Negative", "Splice_Junction_Positive",
                              "Intron_Retention_Negative", "Intron_Retention_Positive"])


    with open(args.output_file + ".mutation_list.overlap.bed", 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            mut_chr, mut_start, mut_end, mut_ref, mut_alt = F[0], int(F[1]) + 1, int(F[2]), F[3], F[4]
            motif_chr, motif_start, motif_end, motif_type, motif_strand, junc_list = F[5], int(F[6]) + 1, int(F[7]) + 1, F[9], F[10], F[11]
            motif_gene = F[8]

            # generate template sequence (splicing junction with and/or without mutation, intron retention with and/or withoug mutation)
            allele_count.generate_template_seq(args.output_file + ".tmp.template_seq.fa" + str(cnum),
                                               args.reference, mut_chr, mut_start, mut_end, mut_ref, mut_alt,
                                               motif_chr, motif_start, motif_end, motif_type, motif_strand,
                                               junc_list, args.donor_size, args.acceptor_size, args.template_size)

            # extract short reads from bam files around the exon-intron boundary
            allele_count.extract_read_around_boundary(args.bam_file, args.output_file + ".tmp.read_seq.fa" + str(cnum),
                                                      args.reference, motif_chr, motif_start, motif_end, args.read_search_margin)

            # perform smith-waterman alignment check
            type2count = pyssw.main2(args.output_file + ".tmp.read_seq.fa" + str(cnum), args.output_file + ".tmp.template_seq.fa" + str(cnum), 
                                     4 * args.template_size - args.template_score_margin)

            print >> hout, '\t'.join([motif_gene, mut_chr, str(mut_start), str(mut_end), mut_ref, mut_alt, 
                             motif_chr, str(motif_start), str(motif_end), motif_type, motif_strand,
                             str(type2count["splice_junction_negative"]), str(type2count["splice_junction_positive"]),
                             str(type2count["intron_retention_negative"]), str(type2count["intron_retention_positive"])])
            
            if not args.debug:
                subprocess.call(["rm", "-rf", args.output_file + ".tmp.template_seq.fa" + str(cnum)])
                subprocess.call(["rm", "-rf", args.output_file + ".tmp.read_seq.fa" + str(cnum)])

            cnum = cnum + 1

    hout.close()

    if not args.debug:
        subprocess.call(["rm", "-rf", args.output_file + ".intron_retention_list.bed"])
        subprocess.call(["rm", "-rf", args.output_file + ".mutation_list.bed"])
        subprocess.call(["rm", "-rf", args.output_file + ".mutation_list.overlap.bed"]) 


