#! /usr/bin/env python

import unittest
import os, tempfile, shutil, filecmp
import intron_retention_utils
from .check_download import *

class TestAlleleCount(unittest.TestCase):

    def setUp(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        check_download("https://storage.googleapis.com/friend1ws_package_data/common/GRCh37.fa", \
                       cur_dir + "/resource/reference_genome/GRCh37.fa")

        check_download("https://storage.googleapis.com/friend1ws_package_data/intron_retention_utils/CCLE-HCC1143-RNA-08.chr21_chr22_20percent.Aligned.sortedByCoord.out.bam", \
                       cur_dir + "/resource/bam/CCLE-HCC1143-RNA-08.chr21_chr22_20percent.Aligned.sortedByCoord.out.bam")

        check_download("https://storage.googleapis.com/friend1ws_package_data/intron_retention_utils/CCLE-HCC1143-RNA-08.chr21_chr22_20percent.Aligned.sortedByCoord.out.bam.bai", \
                       cur_dir + "/resource/bam/CCLE-HCC1143-RNA-08.chr21_chr22_20percent.Aligned.sortedByCoord.out.bam.bai")

        self.parser = intron_retention_utils.parser.create_parser()


    def tearDown(self):
       cur_dir = os.path.dirname(os.path.abspath(__file__))
       shutil.rmtree(cur_dir + "/resource/bam")
 

    def test1(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        input_bam = cur_dir + "/resource/bam/CCLE-HCC1143-RNA-08.chr21_chr22_20percent.Aligned.sortedByCoord.out.bam"
        mutation_file = cur_dir + "/data/mutation/CCLE-HCC1143-DNA-08.genomon_mutation.result.txt"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh37.fa"
        output_file = tmp_dir + "/CCLE-HCC1143-RNA-08.chr21_chr22_20percent.ir_allele_count.result.txt"
        answer_file = cur_dir + "/data/allele_count/CCLE-HCC1143-RNA-08.chr21_chr22_20percent.ir_allele_count.result.txt"

        args = self.parser.parse_args(["allele_count", input_bam, mutation_file, output_file, "--reference", ref_genome])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    unittest.main()


