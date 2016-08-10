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

Simple intron retention count program.
Calculate the number of reads covering each exon-intron boundary and 
putative intron retention reads (that covering enlarged reagion by specified margin size (e.g. -5bp and +5bp from that boundary)).

```
genomon_intron_retention simple_count [-h] 
                                      [-q mapping_qual_thres] 
                                      [--chr_name_list chr_name_list.txt]
                                      [--debug] 
                                      sequence.bam output_file refGene.txt.gz
```

#### result

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



