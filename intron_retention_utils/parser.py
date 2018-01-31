#! /usr/bin/env python

from run import *
import argparse

def create_parser():

    parser = argparse.ArgumentParser(prog = "intron_retention_utils")

    parser.add_argument("--version", action = "version", version = "intron_retention_utils-0.4.0")

    subparsers = parser.add_subparsers()

    ##########
    # simple intron retention count 

    simple_count = subparsers.add_parser("simple_count",
                                         help = "simple intron retention count program")

    simple_count.add_argument("bam_file", metavar = "sequence.bam", default = None, type = str,
                              help = "the path to the bam file")

    simple_count.add_argument("output_file", metavar = "output_file", default = None, type = str, 
                              help = "the path to the output")

    simple_count.add_argument("--grc", default = False, action = 'store_true',
                              help = "convert chromosome names to Genome Reference Consortium nomenclature (default: %(default)s)")

    simple_count.add_argument("--genome_id", choices = ["hg19", "hg38", "mm10"], default = "hg19",
                              help = "the genome id used for selecting UCSC-GRC chromosome name corresponding files (default: %(default)s)")

    simple_count.add_argument("--intron_retention_check_size", metavar = "intron_retention_check_size", default = 10, type = int,
                              help = "exon and intron region size to be covered by putative intron reads (default: %(default)s)")

    simple_count.add_argument("--mapping_qual_thres", metavar = "mapping_qual_thres", default= 20, type=int,
                              help = "threshold for mapping quality for calculating base counts (default: %(default)s)")

    simple_count.add_argument("--keep_improper_pair", action = 'store_true', default = False,
                              help = "keep improper paired reads (activate for single end reads) (default: %(default)s)")

    simple_count.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")

    simple_count.set_defaults(func = simple_count_main)

    ##########
    # allele specific intron retention count

    allele_count = subparsers.add_parser("allele_count",
                                         help = "intron retention count incorporating somatic mutations")

    allele_count.add_argument("bam_file", metavar = "sequence.bam", default = None, type = str,
                              help = "the path to the bam file")

    allele_count.add_argument("mutation_file", metavar = "mutation.txt", default = None, type = str,
                              help = "the path to the mutation list file")

    allele_count.add_argument("output_file", metavar = "output.txt", default = None, type = str,
                              help = "the path to output file")

    allele_count.add_argument("--reference", metavar = "reference.fa", default = None, type = str, required = True,
                              help = "the path to the reference genome file")

    allele_count.add_argument("--grc", default = False, action = 'store_true',
                              help = "convert chromosome names to Genome Reference Consortium nomenclature (default: %(default)s)")

    allele_count.add_argument("--genome_id", choices = ["hg19", "hg38", "mm10"], default = "hg19",
                              help = "the genome id used for selecting UCSC-GRC chromosome name corresponding files (default: %(default)s)")

    allele_count.add_argument("--donor_size", metavar = "donor_size", default = "2,6", type = str,
                              help = "splicing donor site size (exonic region size, intronic region size) (default: %(default)s)")

    allele_count.add_argument("--acceptor_size", metavar = "acceptor_size", default = "8,1", type = str,
                              help = "splicing donor site size (intronic region size, exonic region size) (default: %(default)s)")

    allele_count.add_argument("--template_size", metavar = "check_size", default = 10, type = int,
                              help = "the template sequence sizes for checking intron retention (default: %(default)s)")

    allele_count.add_argument("--template_score_margin", metavar = "check_size", default = 3, type = int,
                              help = "the margin size for checking template match score (default: %(default)s)")

    allele_count.add_argument("--read_search_margin", metavar = "read_search_margin", default = 10, type = int,
                              help = "margin size for extracting short reads around exon-intron junctions (default: %(default)s)")

    allele_count.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")

    allele_count.set_defaults(func = allele_count_main)

    ##########
    # merge control 
    merge_control = subparsers.add_parser("merge_control",
                                          help = "merge, compress and index control intron retention files generated by simple_count command")

    merge_control.add_argument("intron_retention_list", metavar = "intron_retention_list.txt", default = None, type = str,
                               help = "intron retention file path list")

    merge_control.add_argument("output_file", metavar = "merge_control.bedpe.gz", default = None, type = str,
                               help = "the path of the bgzip-compressed and tabix-indexed output file")


    merge_control.add_argument("--ratio_thres", type = float, default = 0.05,
                               help = "register intron retentions whose ratios (Intron_Retention_Read_Count / Edge_Read_Count) \
                               are above this value in at least the specified number of files by sample_num_thres (default: %(default)s)")

    merge_control.add_argument("--sample_num_thres", type = int, default = 2,
                               help = "register intron retentions whose ratios (Intron_Retention_Read_Count / Edge_Read_Count) \
                               are above the specifed value by ratio_thres in at least the specified number of sample by this parameter \
                                (default: %(default)s)")

    merge_control.set_defaults(func = merge_control_main)

    ##########
    # filter

    filter = subparsers.add_parser("filter",
                                   help = "filter out intron retentions that do not satisty specified conditions")

    filter.add_argument("intron_retention_file", metavar = "intron_retention.txt", default = None, type = str,
                        help = "the path to intron retention file generated by simple_count command")

    filter.add_argument("output_file", metavar = "output.txt", default = None, type = str,
                        help = "the path to the output file")

    filter.add_argument("--num_thres", type = int, default = 3,
                        help = "remove intron retentions whose supporting read numbers are below this value (default: %(default)s)")

    filter.add_argument("--ratio_thres", type = int, default = 0.05,
                        help = "remove intron retentions whose whose ratios (Intron_Retention_Read_Count / Edge_Read_Count) \
                        are below this value (default: %(default)s)")
                               
    filter.add_argument("--pooled_control_file", default = None, type = str,
                        help = "the path to control data created by merge_control (default: %(default)s)")

    filter.set_defaults(func = filter_main)

    ##########
    # associate

    associate = subparsers.add_parser("associate",
                                      help = "associate intron retentions with mutations")

    associate.add_argument("intron_retention_file", metavar = "intron_retention.txt", default = None, type = str,
                           help = "the path to intron retention file generated by simple_count command")

    associate.add_argument("mutation_file", metavar = "mutation.txt", default = None, type = str,
                           help = "the path to the mutation list file (vcf or annovar format)")

    associate.add_argument("output_file", metavar = "output_file", default = None, type = str, 
                           help = "the path to the output file")

    associate.add_argument("--donor_size", metavar = "donor_size", default = "3,6", type = str,
                           help = "splicing donor site size (exonic region size, intronic region size) (default: %(default)s)")

    associate.add_argument("--acceptor_size", metavar = "acceptor_size", default = "6,1", type = str,
                           help = "splicing donor site size (intronic region size, exonic region size) (default: %(default)s)")

    associate.add_argument('--mutation_format', choices=['vcf', 'anno'], default = 'vcf',
                           help = "the format of mutation file vcf or annovar (tsv) format (default: %(default)s)")

    associate.add_argument("--reference", metavar = "reference.fa", type = str, 
                           help = "the path to the reference genomoe sequence (necessary when --mutation format is anno)")

    associate.add_argument('--sv', action='store_true',
                           help = "analysis structural variation file")

    associate.add_argument("--intron_margin", metavar = "intron_margin", default = 10, type = int,
                           help = "the margin size of intron region for checking associations with sv breakpoints (default: %(default)s)")

    associate.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")

    associate.set_defaults(func = associate_main)
    ##########

    return parser

