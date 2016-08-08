#! /usr/bin/env python

import sys, os, subprocess
import utils, intron_db, mutation

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
    subprocess.call(["bedtools", "intersect", "-a", args.output_file + ".intron_retention_list.bed",
                     "-b", args.output_file + ".mutation_list.bed", "-wa", "-wb"], stdout = hout)
    hout.close()

