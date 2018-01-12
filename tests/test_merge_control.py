#! /usr/bin/env python

import unittest
import os, glob, tempfile, shutil, filecmp
import intron_retention_utils
from check_download import *

class TestMergeControl(unittest.TestCase):

    def setUp(self):
        self.parser = intron_retention_utils.parser.create_parser()


    def tearDown(self):
        pass 


    def test1(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        all_simple_count_file = glob.glob(cur_dir + "/data/simple_count/*.ir_simple_count.result.txt")
        with open(tmp_dir + "/CCLE.ir_simple_count.reulst_list.txt", 'w') as hout:
            for simple_count_file in sorted(all_simple_count_file):
                print >> hout, simple_count_file

        input_list_file = tmp_dir + "/CCLE.ir_simple_count.reulst_list.txt"
        output_file = tmp_dir + "/merge_control.bed.gz"
        answer_file = cur_dir + "/data/merge_control/merge_control.bed.gz"
 
        args = self.parser.parse_args(["merge_control", input_list_file, output_file])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    unittest.main()


