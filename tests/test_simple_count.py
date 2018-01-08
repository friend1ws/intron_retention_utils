#! /usr/bin/env python

import unittest
import os, tempfile, shutil, filecmp
import intron_retention_utils

class TestSimpleCount(unittest.TestCase):

    def setUp(self):
        self.parser = intron_retention_utils.parser.create_parser()

 
    def test1(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        input_bam = cur_dir + "/data/MCF-7.Aligned.sortedByCoord.out.bam"
        output_file = tmp_dir + "/MCF-7.ir_simple_count.txt"
        answer_file = cur_dir + "/data/simple_count/MCF-7.ir_simple_count.txt"
 
        args = self.parser.parse_args(["simple_count", input_bam, output_file, "--grc"])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    unittest.main()


