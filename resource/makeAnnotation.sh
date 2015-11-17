#! /bin/sh


perl coding_RefSeq.pl refGene.txt refGene.coding.exon.bed refGene.coding.intron.bed refGene.coding.5putr.bed refGene.coding.3putr.bed

mergeBed -i refGene.coding.exon.bed | sortBed -i stdin > refGene.merged.coding.exon.bed

mergeBed -i refGene.coding.intron.bed | sortBed -i stdin > refGene.merged.coding.intron.bed

mergeBed -i refGene.coding.5putr.bed | sortBed -i stdin > refGene.merged.coding.5putr.bed

mergeBed -i refGene.coding.3putr.bed | sortBed -i stdin > refGene.merged.coding.3putr.bed


perl noncoding_RefSeq.pl refGene.txt refGene.noncoding.exon.bed refGene.noncoding.intron.bed

mergeBed -i refGene.noncoding.exon.bed | sortBed -i stdin > refGene.merged.noncoding.exon.bed

mergeBed -i refGene.noncoding.intron.bed | sortBed -i stdin > refGene.merged.noncoding.intron.bed


perl -F"\t" -a -n -e 'print $F[2] . "\t" . $F[4] . "\t" . $F[5] . "\n"' refGene.txt | mergeBed -i stdin | complementBed -i stdin -g ../hg19.genome | sortBed -i stdin > refGene.merged.intergene.bed

perl mergeFiles.pl > refGene.annotation.bed
