# GenomonIntronRetention

## Introduction

GenomonIntronRetention is a software for calculating intron retention events genome-wide from RNA sequencing data.

## Dependency

### Python

Python (>= 2.7), `pysam (>= 0.8.1)` packages

### Software

bedtools (>= 2.20.0)

## Install 
```
git clone  https://github.com/Genomon-Project/GenomonIntronRetention.git
cd GenomonIntronRetention
python setup.py build install
```

## Preparation

For **simple_count** and **allele_count** commands, `refGene.txt.gz` file from UCSC is necessary:
```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
```

Also, if your BAM file is aligned to reference genome using name convention other than UCSC,
you need to set up correspondance file. For example, for making hg19-GRCh37 correspondance file):
```
cd resource
bash make_ucsc_grch.sh
```

Also, for **allele_count** command, a [Smith-Waterman shared library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library) from Mengyao Zhao is necessary.
Create the `libssw.so` and add the path to the LD_LIBRARY_PATH environment variable.



## Commands

### simple_count
```
genomon_intron_retention simple_count [-h] 
                                      [-q mapping_qual_thres] 
                                      [--chr_name_list chr_name_list.txt]
                                      [--debug] 
                                      sequence.bam output_file refGene.txt.gz
```

#### result

1. chromosome name
1. start coordinate
1. end coordinate
1. exon-intron junction name
1. score (currently no meaning)
1. direction of junction
1. the number of read at the junction
1. the number of read identified as intron retention

### allele_count

```
genomon_intron_retention allele_count [-h] 
                                      [--donor_size donor_size]
                                      [--acceptor_size acceptor_size] 
                                      [--chr_name_list chr_name_list.txt] 
                                      [--template_size check_size] 
                                      [--template_score_margin check_size]
                                      [--read_search_margin read_search_margin]
                                      [--debug]
                                      sequence.bam mutation.txt output.txt reference.fa refGene.txt.gz
```
genomon_intron_retention [-h] [--version] [-q mapping_qual_thres] sequence.bam output_prefix annotation_dir bedtools_path
```



