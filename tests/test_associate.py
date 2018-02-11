#! /usr/bin/env python

import unittest
import os, tempfile, shutil, filecmp
import intron_retention_utils
from check_download import *

class TestAssociate(unittest.TestCase):

    def setUp(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        check_download("https://storage.googleapis.com/friend1ws_package_data/common/GRCh37.fa", \
                       cur_dir + "/resource/reference_genome/GRCh37.fa")

        self.parser = intron_retention_utils.parser.create_parser()


    def tearDown(self):
        pass 


    def test1(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        simple_count_file = cur_dir + "/data/simple_count/CCLE-HCC1954-RNA-08.ir_simple_count.result.txt"
        mutation_file = cur_dir + "/data/mutation/CCLE-HCC1954-DNA-08.genomon_mutation.result.txt" 
        output_file = tmp_dir + "/CCLE-HCC1954-RNA-08.ir_simple_count.associate.txt"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh37.fa"
        answer_file = cur_dir + "/data/associate/CCLE-HCC1954-RNA-08.ir_simple_count.associate.txt"
 
        args = self.parser.parse_args(["associate", simple_count_file, mutation_file, output_file, \
                                       "--mutation_format", "anno", "--reference", ref_genome])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)


    def test2(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        simple_count_file = cur_dir + "/data/simple_count/CCLE-GSU-RNA-08.ir_simple_count.txt"
        sv_file = cur_dir + "/data/sv/CCLE-GSU-DNA-08.genomonSV.result.filt.txt"
        output_file = tmp_dir + "/CCLE-GSU-RNA-08.ir_simple_count.associate.txt"
        answer_file = cur_dir + "/data/associate/CCLE-GSU-RNA-08.ir_simple_count.associate.txt"

        args = self.parser.parse_args(["associate", "--sv", simple_count_file, sv_file, output_file])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)


if __name__ == "__main__":
    unittest.main()


