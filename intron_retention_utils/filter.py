#! /usr/bin/env python

from __future__ import print_function
import pysam

def filter_intron_retention(input_file, output_file, control_file, num_thres, ratio_thres):

    is_control = True if control_file is not None else False
    if is_control: control_db = pysam.TabixFile(control_file)

    header2ind = {}
    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:

        header = hin.readline().rstrip('\n').split('\t')
        for i, cname in enumerate(header):
            header2ind[cname] = i

        print('\t'.join(header), file = hout)

        for line in hin:

            F = line.rstrip('\n').split('\t')
            if int(F[header2ind["Intron_Retention_Read_Count"]]) < num_thres: continue

            intron_ratio = 0
            if F[header2ind["Edge_Read_Count"]] != "0":
                intron_ratio = float(F[header2ind["Intron_Retention_Read_Count"]]) / float(F[header2ind["Edge_Read_Count"]])

            if intron_ratio < ratio_thres: continue


            ##########
            # remove control files
            if is_control:
                tabixErrorFlag = 0
                try:
                    records = control_db.fetch(F[0], int(F[1]) - 5, int(F[1]) + 5)
                except Exception as inst:
                    # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                    # tabixErrorMsg = str(inst.args)
                    tabixErrorFlag = 1

                control_flag = 0;
                if tabixErrorFlag == 0:
                    for record_line in records:
                        record = record_line.split('\t')
                        if F[0] == record[0] and F[1] == record[1] and F[2] == record[2]:
                            control_flag = 1

                if control_flag == 1: continue

            print('\t'.join(F), file = hout)
 

    hout.close()


