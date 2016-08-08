#! /bin/sh

# create corresponding table for GRCh name and UCSC name
rm -rf GCF_000001405.13.assembly.txt
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.13.assembly.txt
python make_ucsc_grch.py GCF_000001405.13.assembly.txt ucsc2grch.txt

rm -rf refGene.txt.gz*
rm -rf ensGene.txt.gz*

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/ensGene.txt.gz

echo "python list_edges.py | sort -k1,1 -k2,2n -k3,3n > edge.bed"
python list_edges.py | sort -k1,1 -k2,2n -k3,3n > edge.bed

echo "python broaden_junc.py edge.bed 5 | sort -k1,1 -k2,2n -k3,3n > edge_broaden.bed"
python broaden_junc.py edge.bed 5 | sort -k1,1 -k2,2n -k3,3n > edge_broaden.bed

