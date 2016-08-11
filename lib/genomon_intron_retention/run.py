#! /usr/bin/env python

import sys, os, subprocess
import intron_db

def simple_count_main(args):

    import simple_count

    output_dir = os.path.dirname(args.output_file)
    if output_dir != "" and not os.path.exists(output_dir):
       os.makedirs(output_dir)


    # intron_db.generate_edge_bed(args.ref_gene_file, output_prefix + ".refGene.edge.bed", args.chr_name_list)
    intron_db.generate_intron_retention_list(args.ref_gene_file, args.output_file + ".refGene.edge.bed", 
                                             "1,0", "0,1", args.chr_name_list)

    intron_db.broaden_edge(args.output_file + ".refGene.edge.bed", args.output_file + ".refGene.edge_broaden.bed", 5)

    simple_count.filterImproper(args.bam_file, args.output_file + ".filt.bam", args.q)


    hout = open(args.output_file + ".filt.bed12", 'w')
    s_ret = subprocess.call(["bedtools", "bamtobed", "-bed12", "-i", args.output_file + ".filt.bam"], stdout = hout)
    hout.close()

    if s_ret != 0:
        print >> sys.stderr, "error in generating filt.bed12"
        sys.exit(1)


    hout = open(args.output_file + ".edge.bed", 'w')
    s_ret = subprocess.call(["bedtools", "intersect", "-a", args.output_file + ".filt.bed12", "-b", 
                             args.output_file + ".refGene.edge.bed", "-split", "-wo"], stdout = hout) 
    hout.close()

    if s_ret != 0:
        print >> sys.stderr, "error in generating edge.bed"
        sys.exit(1)


    hout = open(args.output_file + ".edge_broaden.bed", 'w')
    s_ret = subprocess.call(["bedtools", "intersect", "-a", args.output_file + ".filt.bed12", "-b", 
                             args.output_file + ".refGene.edge_broaden.bed", "-split", "-wo"], stdout = hout) 
    hout.close()

    if s_ret != 0:
        print >> sys.stderr, "error in generating edge_broaden.bed"
        sys.exit(1)

    simple_count.summarize_edge(args.output_file + ".edge.bed", args.output_file + ".edge_broaden.bed", args.output_file + ".unsorted", 5)

    # print header
    hout = open(args.output_file, 'w')
    print >> hout, '\t'.join(["Chr", "Boundary_Pos", "Gene_Symbol", "Motif_Type", "Strand", 
                              "Junction_List", "Gene_ID_List", "Exon_Num_List", "Edge_Read_Count", "Intron_Retention_Read_Count"])
    hout.close()

    # sort the result
    hout = open(args.output_file, 'a')
    subprocess.call(["sort", "-k1,1", "-k2,2n", "-k3,3n", args.output_file + ".unsorted"], stdout = hout)
    hout.close()

    if not args.debug:
        subprocess.call(["rm", "-rf", args.output_file + ".filt.bam"])
        subprocess.call(["rm", "-rf", args.output_file + ".filt.bed12"])
        subprocess.call(["rm", "-rf", args.output_file + ".refGene.edge.bed"])
        subprocess.call(["rm", "-rf", args.output_file + ".refGene.edge_broaden.bed"])
        subprocess.call(["rm", "-rf", args.output_file + ".edge.bed"])
        subprocess.call(["rm", "-rf", args.output_file + ".edge_broaden.bed"])
        subprocess.call(["rm", "-rf", args.output_file + ".unsorted"])


def allele_count_main(args):

    import mutation, allele_count, pyssw

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


def merge_control_main(args):

    # make directory for output if necessary
    if os.path.dirname(args.output_file) != "" and not os.path.exists(os.path.dirname(args.output_file)):
        os.makedirs(os.path.dirname(args.output_file))

    hin = open(args.intron_retention_list, 'r')
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
       
                    key = '\t'.join(F[:8]) 
                    if intron_ratio >= args.ratio_thres: 
                        print >> hout, key + '\t' + str(round(intron_ratio, 3)) + '\t' + read_count
                

    hout = open(args.output_file + ".sorted", 'w')
    s_ret = subprocess.call(["sort", "-k1,1", "-k2,2n", args.output_file + ".unsorted"], stdout = hout)
    hout.close()

    if s_ret != 0:
        print >> sys.stderr, "Error in sorting merged junction file"
        sys.exit(1)

    hout = open(args.output_file + ".merged", 'w')
    with open(args.output_file + ".sorted", 'r') as hin:
        temp_key = ""
        temp_ratio = []
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = '\t'.join(F[:8])
            ratio = F[8]
            read_count = F[9]
            if key != temp_key:
                if temp_key != "":
                    if len(temp_ratio) >= args.sample_num_thres:
                        print >> hout, temp_key + '\t' + ';'.join(temp_ratio) + '\t' + ';'.join(temp_read_count) 
                temp_key = key
                temp_ratio = []
                temp_read_count = []
            else:
                temp_ratio.append(str(ratio))
                temp_read_count.append(read_count)

        if key != temp_key:
            if temp_key != "":
                if len(temp_ratio) >= sample_num_thres:
                    print >> hout, temp_key + '\t' + ';'.join(temp_ratio) + '\t' + ';'.join(temp_read_count)


    hout = open(args.output_file, 'w')
    s_ret = subprocess.call(["bgzip", "-f", "-c", args.output_file + ".merged"], stdout = hout)
    hout.close()

    if s_ret != 0:
        print >> sys.stderr, "Error in compression merged junction file"
        sys.exit(1)


    s_ret = subprocess.call(["tabix", "-p", "vcf", args.output_file])
    if s_ret != 0:
        print >> sys.stderr, "Error in indexing merged junction file"
        sys.exit(1)

    """
    subprocess.call(["rm", "-f", args.output_file + ".unsorted"])
    subprocess.call(["rm", "-f", args.output_file + ".sorted"])
    subprocess.call(["rm", "-f", args.output_file + ".merged"])
    """


def filter_main(args):

    import filter
    filter.filter_intron_retention(args.intron_retention_file, args.output_file, args.pooled_control_file, args.num_thres, args.ratio_thres)


def associate_main(args):
    
    import associate
    associate.generate_mutation_target(args.intron_retention_file, args.output_file + ".target_list.bed",
                                       args.output_file + ".intron_retention_file.header", args.donor_size, args.acceptor_size)

