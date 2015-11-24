#! /bin/sh

# rm -rf refGene.txt.gz*
# rm -rf ensGene.txt.gz*

# wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
# wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/ensGene.txt.gz

# echo "python list_edges.py | sort -k1,1 -k2,2 -k3,3 > edge.bed"
# python list_edges.py | sort -k1,1 -k2,2 -k3,3 > edge.bed

echo "python broaden_junc.py edge.bed 5 | sort -k1,1 -k2,2 -k3,3 > edge_broaden.bed"
python broaden_junc.py edge.bed 5 | sort -k1,1 -k2,2 -k3,3 > edge_broaden.bed

