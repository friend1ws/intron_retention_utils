#! /usr/bin/env python

from __future__ import print_function
import sys, os, subprocess

import annot_utils.boundary

from . import intron_db
from .logger import get_logger
logger = get_logger()

try:
    from builtins import round
except ImportError:
    pass


def simple_count_main(args):

    from . import simple_count
    from annot_utils.utils import grc_check

    if args.grc == True:
        logger.warning("--grc argument is deprecated and ignored.")

    is_grc = grc_check(args.bam_file)
    
    output_dir = os.path.dirname(args.output_file)
    if output_dir != "" and not os.path.exists(output_dir):
       os.makedirs(output_dir)


    # intron_db.generate_edge_bed(args.ref_gene_file, output_prefix + ".refGene.edge.bed", args.chr_name_list)
    # intron_db.generate_intron_retention_list(args.ref_gene_file, args.output_file + ".refGene.edge.bed", 
    #                                          "1,0", "0,1", args.chr_name_list)

    annot_utils.boundary.make_boundary_info(args.output_file + ".refGene.edge.bed.gz", args.genome_id, is_grc, "1,0", "0,1")

    intron_db.broaden_edge(args.output_file + ".refGene.edge.bed.gz", 
                           args.output_file + ".refGene.edge_broaden.bed",
                           args.intron_retention_check_size)

    simple_count.filterImproper(args.bam_file, args.output_file + ".filt.bam", args.mapping_qual_thres, args.keep_improper_pair)


    hout = open(args.output_file + ".filt.bed12", 'w')
    s_ret = subprocess.check_call(["bedtools", "bamtobed", "-bed12", "-i", args.output_file + ".filt.bam"], stdout = hout)
    hout.close()

    if s_ret != 0:
        logger.error("Error in generating filt.bed12.")
        sys.exit(1)


    hout = open(args.output_file + ".edge.bed", 'w')
    s_ret = subprocess.check_call(["bedtools", "intersect", "-a", args.output_file + ".filt.bed12", "-b", 
                             args.output_file + ".refGene.edge.bed.gz", "-split", "-wo"], stdout = hout) 
    hout.close()

    if s_ret != 0:
        logger.error("Error in generating edge.bed.")
        sys.exit(1)


    hout = open(args.output_file + ".edge_broaden.bed", 'w')
    s_ret = subprocess.check_call(["bedtools", "intersect", "-a", args.output_file + ".filt.bed12", "-b", 
                             args.output_file + ".refGene.edge_broaden.bed", "-split", "-wo"], stdout = hout) 
    hout.close()

    if s_ret != 0:
        logger.error("Error in generating edge_broaden.bed")
        sys.exit(1)

    simple_count.summarize_edge(args.output_file + ".edge.bed",
                                args.output_file + ".edge_broaden.bed",
                                args.output_file + ".unsorted",
                                args.intron_retention_check_size)

    # print header
    hout = open(args.output_file, 'w')
    print('\t'.join(["Chr", "Boundary_Pos", "Gene_Symbol", "Motif_Type", "Strand", 
        "Junction_List", "Gene_ID_List", "Exon_Num_List", "Edge_Read_Count", "Intron_Retention_Read_Count"]), file = hout)
    hout.close()

    # sort the result
    hout = open(args.output_file, 'a')
    subprocess.check_call(["sort", "-f", "-k1,1", "-k2,2n", "-k3,3", args.output_file + ".unsorted"], stdout = hout)
    hout.close()

    if not args.debug:
        subprocess.check_call(["rm", "-rf", args.output_file + ".filt.bam"])
        subprocess.check_call(["rm", "-rf", args.output_file + ".filt.bed12"])
        subprocess.check_call(["rm", "-rf", args.output_file + ".refGene.edge.bed.gz"])
        subprocess.check_call(["rm", "-rf", args.output_file + ".refGene.edge.bed.gz.tbi"])
        subprocess.check_call(["rm", "-rf", args.output_file + ".refGene.edge_broaden.bed"])
        subprocess.check_call(["rm", "-rf", args.output_file + ".edge.bed"])
        subprocess.check_call(["rm", "-rf", args.output_file + ".edge_broaden.bed"])
        subprocess.check_call(["rm", "-rf", args.output_file + ".unsorted"])


def allele_count_main(args):

    from . import mutation, allele_count, pyssw

    from annot_utils.utils import grc_check

    if args.grc == True:
        logger.warning("--grc argument is deprecated and ignored.")

    is_grc = grc_check(args.bam_file)


    output_dir = os.path.dirname(args.output_file)
    if output_dir != "" and not os.path.exists(output_dir):
       os.makedirs(output_dir)

    # intron_db.generate_intron_retention_list(args.ref_gene_file, args.output_file + ".intron_retention_list.bed", 
    #                                          args.donor_size, args.acceptor_size, args.chr_name_list)

    annot_utils.boundary.make_boundary_info(args.output_file + ".refGene.edge.bed.gz", args.genome_id, is_grc, args.donor_size, args.acceptor_size)

    mutation.anno2bed(args.mutation_file, args.output_file + ".mutation_list.bed")

    hout = open(args.output_file + ".mutation_list.overlap.bed", 'w')
    subprocess.check_call(["bedtools", "intersect", "-a", args.output_file + ".mutation_list.bed",
                     "-b", args.output_file + ".refGene.edge.bed.gz", "-wa", "-wb"], stdout = hout)
    hout.close()
   
    """ 
    cnum = 0
    hout = open(args.output_file, 'w')
    # print header
    print('\t'.join(["Gene_Symbol", "Chr_Mut", "Start_Mut", "End_Mut", "Ref_Mut", "Alt_Mut", 
                     "Chr_Motif", "Start_Motif", "End_Motif", "Type_Motif", "Strand_Motif",
                     "Splice_Junction_Negative", "Splice_Junction_Positive",
                     "Intron_Retention_Negative", "Intron_Retention_Positive"]), file = hout)
    """

    cnum = 0
    hout = open(args.output_file + ".unsorted", 'w')
    with open(args.output_file + ".mutation_list.overlap.bed", 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            mut_chr, mut_start, mut_end, mut_ref, mut_alt = F[0], int(F[1]) + 1, int(F[2]), F[3], F[4]
            motif_chr, motif_start, motif_end, motif_type, motif_strand, junc_list = F[5], int(F[6]) + 1, int(F[7]), F[9], F[10], F[11]
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

            print('\t'.join([motif_gene, mut_chr, str(mut_start), str(mut_end), mut_ref, mut_alt, 
                motif_chr, str(motif_start), str(motif_end), motif_type, motif_strand,
                str(type2count["splice_junction_negative"]), str(type2count["splice_junction_positive"]),
                str(type2count["intron_retention_negative"]), str(type2count["intron_retention_positive"])]), file = hout)
            
            if not args.debug:
                subprocess.check_call(["rm", "-rf", args.output_file + ".tmp.template_seq.fa" + str(cnum)])
                subprocess.check_call(["rm", "-rf", args.output_file + ".tmp.read_seq.fa" + str(cnum)])

            cnum = cnum + 1

    hout.close()


    hout = open(args.output_file, 'w')
    # print header
    print('\t'.join(["Gene_Symbol", "Chr_Mut", "Start_Mut", "End_Mut", "Ref_Mut", "Alt_Mut", 
                     "Chr_Motif", "Start_Motif", "End_Motif", "Type_Motif", "Strand_Motif",
                     "Splice_Junction_Negative", "Splice_Junction_Positive",
                     "Intron_Retention_Negative", "Intron_Retention_Positive"]), file = hout)
    hout.close()


    # sort the result
    hout = open(args.output_file, 'a')
    subprocess.check_call(["sort", "-f", "-k2,2", "-k3,3n", "-k1,1", args.output_file + ".unsorted"], stdout = hout)
    hout.close()


    if not args.debug:
        subprocess.check_call(["rm", "-rf", args.output_file + ".refGene.edge.bed.gz"])
        subprocess.check_call(["rm", "-rf", args.output_file + ".refGene.edge.bed.gz.tbi"])
        subprocess.check_call(["rm", "-rf", args.output_file + ".mutation_list.bed"])
        subprocess.check_call(["rm", "-rf", args.output_file + ".mutation_list.overlap.bed"]) 
        subprocess.check_call(["rm", "-rf", args.output_file + ".unsorted"])

def merge_control_main(args):

    # make directory for output if necessary
    if os.path.dirname(args.output_file) != "" and not os.path.exists(os.path.dirname(args.output_file)):
        os.makedirs(os.path.dirname(args.output_file))

    hout = open(args.output_file + ".unsorted", 'w')
    header2ind = {}
    with open(args.intron_retention_list, 'r') as hin:
        for line in hin:

            intron_retention_file = line.rstrip('\n')
            with open(intron_retention_file, 'r') as hin2:
           
                header = hin2.readline().rstrip('\n').split('\t')
                for i, cname in enumerate(header):
                    header2ind[cname] = i

                for line2 in hin2:

                    F = line2.rstrip('\n').split('\t')
                    intron_ratio = 0
                    read_count = F[header2ind["Edge_Read_Count"]] + ',' + F[header2ind["Intron_Retention_Read_Count"]]
                    if F[header2ind["Edge_Read_Count"]] != "0":
                        intron_ratio = float(F[header2ind["Intron_Retention_Read_Count"]]) / float(F[header2ind["Edge_Read_Count"]])
       
                    key = intron_retention_file + '\t' + '\t'.join(F[:8]) 
                    if intron_ratio >= args.ratio_thres: 
                        print(key + '\t' + str(round(intron_ratio, 3)) + '\t' + read_count, file = hout)
    hout.close()            

    hout = open(args.output_file + ".sorted", 'w')
    s_ret = subprocess.check_call(["sort", "-f", "-k2,2", "-k3,3n", "-k4,4", "-k1,1", args.output_file + ".unsorted"], stdout = hout)
    hout.close()

    if s_ret != 0:
        logger.error("Error in sorting merged junction file.")
        sys.exit(1)

    hout = open(args.output_file + ".merged", 'w')
    with open(args.output_file + ".sorted", 'r') as hin:
        temp_key = ""
        temp_ratio = []
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = '\t'.join(F[1:9])
            ratio = F[9]
            read_count = F[10]
            if key != temp_key:
                if temp_key != "":
                    if len(temp_ratio) >= args.sample_num_thres:
                        print(temp_key + '\t' + ';'.join(temp_ratio) + '\t' + ';'.join(temp_read_count), file = hout)

                temp_key = key
                temp_ratio = []
                temp_read_count = []
                
            temp_ratio.append(str(ratio))
            temp_read_count.append(read_count)

        if key != temp_key:
            if temp_key != "":
                if len(temp_ratio) >= sample_num_thres:
                    print(temp_key + '\t' + ';'.join(temp_ratio) + '\t' + ';'.join(temp_read_count), file = hout)
    hout.close()


    hout = open(args.output_file, 'w')
    s_ret = subprocess.check_call(["bgzip", "-f", "-c", args.output_file + ".merged"], stdout = hout)
    hout.close()

    if s_ret != 0:
        logger.error("Error in compression merged intron retention file.")
        sys.exit(1)


    s_ret = subprocess.check_call(["tabix", "-p", "vcf", args.output_file])
    if s_ret != 0:
        logger.error("Error in indexing merged intron retention file.")
        sys.exit(1)

    subprocess.check_call(["rm", "-f", args.output_file + ".unsorted"])
    subprocess.check_call(["rm", "-f", args.output_file + ".sorted"])
    subprocess.check_call(["rm", "-f", args.output_file + ".merged"])


def filter_main(args):

    from . import filter
    filter.filter_intron_retention(args.intron_retention_file, args.output_file, args.pooled_control_file, args.num_thres, args.ratio_thres)


def associate_main(args):
    
    from . import mutation, associate

    is_sv = True if args.sv else False

    if args.mutation_format == "anno":
        logger.warning("--mutation_format is deprepaced and ignored.")

    is_anno = False if args.mutation_format.endswith(".vcf") or args.mutation_format.endswith(".vcf.gz") else True

    if is_sv == False:

        associate.generate_mutation_target(args.intron_retention_file, args.output_file + ".target_list.bed",
                                           args.output_file + ".intron_retention_file.header", args.donor_size, args.acceptor_size)

        # convert annovar format to vcf
        if is_anno:
            mutation.anno2vcf(args.mutation_file, args.output_file + ".tmp.mutation.unsorted.vcf", args.reference)
        else:
            mutation.remove_vcf_header(args.mutation_file, args.output_file + ".tmp.mutation.unsorted.vcf")

        hout = open(args.output_file + ".tmp.mutation.sorted.vcf", 'w')
        s_ret = subprocess.check_call(["sort", "-f", "-k1,1", "-k2,2n", args.output_file + ".tmp.mutation.unsorted.vcf"], stdout = hout)
        hout.close()

        if s_ret != 0:
            logger.error("Error in sorting vcf file.")
            sys.exit(1)

        mutation.vcf2bed(args.output_file + ".tmp.mutation.sorted.vcf", args.output_file + ".tmp.mutation.sorted.vcf.bed")

        # mutation.anno2bed(args.mutation_file, args.output_file + ".mutation_list.bed")
       
        hout = open(args.output_file + ".mutation_list.associate.bed", 'w')
        subprocess.check_call(["bedtools", "intersect", "-a", args.output_file + ".tmp.mutation.sorted.vcf.bed",
                         "-b", args.output_file + ".target_list.bed", "-wa", "-wb"], stdout = hout)
        hout.close()

        associate.process_result(args.output_file + ".mutation_list.associate.bed", 
                                 args.output_file + ".intron_retention_file.header", 
                                 args.output_file, args.donor_size, args.acceptor_size)

        if args.debug == False:
            subprocess.check_call(["rm", "-rf", args.output_file + ".target_list.bed"])
            subprocess.check_call(["rm", "-rf", args.output_file + ".intron_retention_file.header"])
            subprocess.check_call(["rm", "-rf", args.output_file + ".tmp.mutation.unsorted.vcf"])
            subprocess.check_call(["rm", "-rf", args.output_file + ".tmp.mutation.sorted.vcf"])
            subprocess.check_call(["rm", "-rf", args.output_file + ".tmp.mutation.sorted.vcf.bed"])
            subprocess.check_call(["rm", "-rf", args.output_file + ".mutation_list.associate.bed"])


    else:
        
        associate.generate_sv_target(args.intron_retention_file, args.output_file + ".target_list.bed",
                                     args.output_file + ".intron_retention_file.header", args.intron_margin)

        mutation.genosv2bed(args.mutation_file, args.output_file + ".tmp.unsorted.bedpe")

        hout = open(args.output_file + ".tmp.bedpe", 'w')
        s_ret = subprocess.check_call(["sort", "-k1,1", "-k2,2n", "-k3,3n", "-k4,4", "-k5,5n", "-k6,6n", args.output_file + ".tmp.unsorted.bedpe"], stdout = hout)
        hout.close()

        if s_ret != 0:
            logger.error("Error in sorting bedpe file.")
            sys.exit(1)

        hout = open(args.output_file + ".sv_list.associate.bed", 'w')
        subprocess.check_call(["bedtools", "intersect", "-a", args.output_file + ".tmp.bedpe",
                         "-b", args.output_file + ".target_list.bed", "-wa", "-wb"], stdout = hout)
        hout.close()

        associate.process_result_sv(args.output_file + ".sv_list.associate.bed",
                                    args.output_file + ".intron_retention_file.header",
                                    args.output_file)

        if args.debug == False:
            subprocess.check_call(["rm", "-rf", args.output_file + ".target_list.bed"])
            subprocess.check_call(["rm", "-rf", args.output_file + ".intron_retention_file.header"])
            subprocess.check_call(["rm", "-rf", args.output_file + ".tmp.unsorted.bedpe"])
            subprocess.check_call(["rm", "-rf", args.output_file + ".tmp.bedpe"])
            subprocess.check_call(["rm", "-rf", args.output_file + ".sv_list.associate.bed"])


