# intron_retention_utils

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://travis-ci.org/friend1ws/intron_retention_utils.svg?branch=devel)](https://travis-ci.org/friend1ws/intron_retention_utils)

A software for calculating intron retention events genome-wide from RNA sequencing data.

## Dependency

### Python

Python (>= 2.7), `pysam`,[`annot_utils`](https://github.com/friend1ws/annot_utils) packages.

### Software

[bedtools](http://bedtools.readthedocs.io/en/latest/), [hstlib](http://www.htslib.org)

## Install 
```
git clone  https://github.com/friend1ws/intron_retention_util.git
cd intron_retention_utils
python setup.py build install
```

## Preparation

For **allele_count** command, a [Smith-Waterman shared library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library) from Mengyao Zhao is necessary.
Also, if your BAM file is aligned to reference genome using name convention other than UCSC,
Create the `libssw.so` and add the path to the LD_LIBRARY_PATH environment variable.

## Commands

### simple_count

Simple intron retention count program.
Calculate the number of reads covering each exon-intron boundary and 
putative intron retention reads (that covering enlarged reagion by specified margin size (e.g. -5bp and +5bp from that boundary).

```
intron_retention_utils simple_count [-h] 
                                    [--grc]
                                    [--genome_id {hg19,hg38,mm10}]
                                    [--intron_retention_check_size intron_retention_check_size]
                                    [--mapping_qual_thres mapping_qual_thres]
                                    [--keep_improper_pair] [--debug]
                                    sequence.bam output_file
```

#### About result

* **Chr**: chromosome of the exon-intron boundary
* **Boundary_Pos**: coordinate of the exon-intron boundary (the last exonic base)
* **Gene_Symbol**: gene symbol from refGene.txt.gz
* **Motif_Type**: splicing donor or acceptor
* **Strand**: transcription starnd of the gene
* **Junction_List**: cannonical splicing junction list from that exon-intron boundary
* **Gene_ID_List**: refGene ID list with that exon-intron boundary
* **Exon_Num_List**: exon numbers for each refGene IDs
* **Edge_Read_Count**: the number of reads covering each exon-intron boundary
* **Intron_Retention_Read_Count**: the number of putative intron retention reads



### allele_count

```
intron_retention_utils allele_count [-h] 
                                    [--grc]
                                    [--genome_id {hg19,hg38,mm10}]
                                    [--donor_size donor_size]
                                    [--acceptor_size acceptor_size]
                                    [--template_size check_size]
                                    [--template_score_margin check_size]
                                    [--read_search_margin read_search_margin]
                                    [--debug]
                                    sequence.bam mutation.txt
                                    output.txt reference.fa
```

#### About result

* **Gene_Symbol**: gene symbol
* **Chr_Mut**: chromosome of the mutation
* **Start_Mut**: start coordinate of the mutation
* **End_Mut**: end coordinate of the the mutation
* **Ref_Mut**: reference allele of the mutation
* **Alt_Mut**: alternative allele of the mutation
* **Chr_Motif**: chromosome of the splicing motif
* **Start_Motif**: start coordinate of the splicing motif
* **End_Motif**: end coordinate of the splicing motif
* **Type_Motif**: donor or acceptor
* **Strand_Motif**: transcription strand of the gene 
* **Splice_Junction_Negative**: the number of normaly spliced reads without the alternative allele
* **Splice_Junction_Positive**: the number of normaly spliced reads with the alternative allele
* **Intron_Retention_Negative**: the number of putative intron retention reads without the alternative allele
* **Intron_Retention_Positive**: the number of putative intron retention reads with the alternative allele

### merge_control

Merge the intron retention file of control data (typically) for later filtering.

```
intron_retention_utils merge_control [-h] 
                                     [--ratio_thres RATIO_THRES]
                                     [--sample_num_thres SAMPLE_NUM_THRES]
                                     intron_retention_list.txt
                                     output_file
```

### filter

Filter out intron retentions that do not satisty specified conditions
```
intron_retention_utils filter [-h] 
                              [--num_thres NUM_THRES]
                              [--ratio_thres RATIO_THRES]
                              [--pooled_control_file POOLED_CONTROL_FILE]
                              intron_retention.txt output.txt
```

### associate

Associate intron retention counts (typically output of simple_count commands) with mutations
```
intron_retention_utils associate [-h] [--donor_size donor_size]
                                        [--acceptor_size acceptor_size]
                                        [--mutation_format {vcf,anno}]
                                        [--reference reference.fa] [--sv]
                                        [--intron_margin intron_margin]
                                        [--debug]
                                        intron_retention.txt mutation.txt
                                        output_file
```

#### About result
The following columns are added to the input files:

* **Mutation_Key**: vcf format mutation aggregated by commas
* **Motif_Pos**: coordinate of motif positions
* **Mutation_Type**: `splicing donor disruption` or `splicing acceptor disruption`
* **Is_Canonical**: whether the mutation is disrupting cannonical splicing motifs (GT-AG) or not
* **Intron_Retention_Type**: `direct-impact` or `opposite-side-impact`
