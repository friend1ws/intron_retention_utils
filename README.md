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

For *simple_count* and *allele_count* commands, `refGene.txt.gz` file from UCSC is necessary:
```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
```

Also, if your BAM file is aligned to reference genome using name convention other than UCSC,
you need to set up correspondance file. For example, for making hg19-GRCh37 correspondance file):
```
cd resource
bash make_ucsc_grch.sh
```

Also, for *allele_count* command, smith waterman shared library is necessary


First, you need to set up the exon-intron junction information.
The easiest way is to use the prepared script.

```
cd resource
bash prepareIntron.sh 
```

## Commands

```
genomon_intron_retention [-h] [--version] [-q mapping_qual_thres] sequence.bam output_prefix annotation_dir bedtools_path
```

## Results

1. chromosome name
1. start coordinate
1. end coordinate
1. exon-intron junction name
1. score (currently no meaning)
1. direction of junction
1. the number of read at the junction
1. the number of read identified as intron retention

