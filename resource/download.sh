#! /bin/sh

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/ensGene.txt.gz

gunzip knownGene.txt.gz
gunzip refGene.txt.gz
gunzip ensGene.txt.gz


