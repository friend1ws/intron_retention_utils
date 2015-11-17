#! /bin/sh

echo "sh download.sh"
sh download.sh

echo "perl listEdges.pl > edge.bed"
perl listEdges.pl > edge.bed

echo "perl brodenJunc.pl edge.bed > edge_broden.bed"
perl brodenJunc.pl edge.bed > edge_broden.bed


# rm -rf exon.bed
# rm -rf exon.proc.tmp.fasta

