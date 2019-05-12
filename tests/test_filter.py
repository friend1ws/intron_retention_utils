#! /usr/bin/env python

import unittest
import os, tempfile, shutil, filecmp
import intron_retention_utils
from .check_download import *

class TestFilter(unittest.TestCase):

    def setUp(self):
        self.parser = intron_retention_utils.parser.create_parser()


    def tearDown(self):
        pass 


    def test1(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        input_file = cur_dir + "/data/simple_count/CCLE-HCC1954-RNA-08.chr21_chr22_20percent.ir_simple_count.result.txt"
        output_file = tmp_dir + "/CCLE-HCC1954-RNA-08.chr21_chr22_20percent.ir_simple_count.filter.result.txt"
        answer_file = cur_dir + "/data/filter/CCLE-HCC1954-RNA-08.chr21_chr22_20percent.ir_simple_count.filter.result.txt"
 
        args = self.parser.parse_args(["filter", input_file, output_file])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    unittest.main()


