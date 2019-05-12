#! /usr/bin/env python

import unittest
import os, tempfile, shutil, filecmp
import intron_retention_utils
from check_download import *

class TestSimpleCount(unittest.TestCase):

    def setUp(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        check_download("https://storage.googleapis.com/friend1ws_package_data/intron_retention_utils/CCLE-HCC1954-RNA-08.chr21_chr22_20percent.Aligned.sortedByCoord.out.bam", \
                       cur_dir + "/resource/bam/CCLE-HCC1954-RNA-08.chr21_chr22_20percent.Aligned.sortedByCoord.out.bam")

        check_download("https://storage.googleapis.com/friend1ws_package_data/intron_retention_utils/CCLE-HCC1954-RNA-08.chr21_chr22_20percent.Aligned.sortedByCoord.out.bam.bai", \
                       cur_dir + "/resource/bam/CCLE-HCC1954-RNA-08.chr21_chr22_20percent.Aligned.sortedByCoord.out.bam.bai")

        self.parser = intron_retention_utils.parser.create_parser()


    def tearDown(self):
       cur_dir = os.path.dirname(os.path.abspath(__file__))
       shutil.rmtree(cur_dir + "/resource/bam")
 

    def test1(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        input_bam = cur_dir + "/resource/bam/CCLE-HCC1954-RNA-08.chr21_chr22_20percent.Aligned.sortedByCoord.out.bam"
        output_file = tmp_dir + "/CCLE-HCC1954-RNA-08.chr21_chr22_20percent.ir_simple_count.result.txt"
        answer_file = cur_dir + "/data/simple_count/CCLE-HCC1954-RNA-08.chr21_chr22_20percent.ir_simple_count.result.txt"

        print(output_file)
        print(answer_file)
 
        args = self.parser.parse_args(["simple_count", input_bam, output_file])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        # shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    unittest.main()


