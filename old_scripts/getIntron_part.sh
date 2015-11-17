#! /bin/sh
#$ -S /bin/sh
#$ -cwd

source ../lib/config.sh
source ../lib/utility.sh
PATH=${BLAT_PATH}:$PATH
referece=${BLAT_REF}

PATH=${BEDTOOLS_PATH}:${PATH}
PATH=${SAMTOOLS_PATH}:${PATH}

SAMPLE=$1
NUM=${SGE_TASK_ID}
# NUM=$2 

SEQDIR=${OUTPUTDIR}/${SAMPLE}/bam
TARGETDIR=${OUTPUTDIR}/${SAMPLE}/intron

INTERVAL=${DBDIR}/interval_list_hg19_nongap

# sleep 
# sh sleep.sh
 

REGION_A=`head -n1 ${INTERVAL}/${NUM}.interval_list | awk '{split($0, ARRAY, "-"); print ARRAY[1]}'`
REGION_B=`tail -n1 ${INTERVAL}/${NUM}.interval_list | awk '{split($0, ARRAY, "-"); print ARRAY[2]}'`
REGION="${REGION_A}-${REGION_B}"
echo ${REGION} 

##########
# : <<'_COMMENT_OUT_'

echo "samtools view -F 1024 -h ${SEQDIR}/sequence.bam ${REGION} >  ${TARGETDIR}/tmp/temp${NUM}.sam"
samtools view -F 1024 -h ${SEQDIR}/sequence.bam ${REGION} > ${TARGETDIR}/tmp/temp${NUM}.sam
check_error $?


echo "perl filterRead.pl ${TARGETDIR}/tmp/temp${NUM}.sam 60 > ${TARGETDIR}/tmp/temp${NUM}.sam.filt"
perl filterRead.pl ${TARGETDIR}/tmp/temp${NUM}.sam 60 > ${TARGETDIR}/tmp/temp${NUM}.sam.filt
check_error $?

echo "samtools view -h -bS ${TARGETDIR}/tmp/temp${NUM}.sam.filt > ${TARGETDIR}/tmp/temp${NUM}.bam"
samtools view -h -bS ${TARGETDIR}/tmp/temp${NUM}.sam.filt > ${TARGETDIR}/tmp/temp${NUM}.bam
check_error $?


echo "bamToBed -i ${TARGETDIR}/tmp/temp${NUM}.bam -bed12 > ${TARGETDIR}/tmp/temp${NUM}.bed12"
bamToBed -i ${TARGETDIR}/tmp/temp${NUM}.bam -bed12 > ${TARGETDIR}/tmp/temp${NUM}.bed12

rm -rf ${TARGETDIR}/tmp/temp${NUM}.sam
rm -rf ${TARGETDIR}/tmp/temp${NUM}.sam.filt 
rm -rf ${TARGETDIR}/tmp/temp${NUM}.bam

# _COMMENT_OUT_
##########

# :<<_COMMENT_OUT_

echo "intersectBed -a ${TARGETDIR}/tmp/temp${NUM}.bed12 -b ${DBDIR}/intron/edge.bed -split -wao > ${TARGETDIR}/tmp/temp${NUM}.edge.bed"
intersectBed -a ${TARGETDIR}/tmp/temp${NUM}.bed12 -b ${DBDIR}/intron/edge.bed -split -wao > ${TARGETDIR}/tmp/temp${NUM}.edge.bed
check_error $?

perl -F"\t" -ane 'if ($F[$#F] > 0) {print join("\t", @F);}' ${TARGETDIR}/tmp/temp${NUM}.edge.bed  | sort -k 4 - > ${TARGETDIR}/tmp/temp${NUM}.edge.sorted.bed
check_error $?

echo "intersectBed -a ${TARGETDIR}/tmp/temp${NUM}.bed12 -b ${DBDIR}/intron/edge_broden.bed -split -wao > ${TARGETDIR}/tmp/temp${NUM}.edge_broden.bed"
intersectBed -a ${TARGETDIR}/tmp/temp${NUM}.bed12 -b ${DBDIR}/intron/edge_broden.bed -split -wao > ${TARGETDIR}/tmp/temp${NUM}.edge_broden.bed
check_error $?

perl -F"\t" -ane 'if ($F[$#F] > 0) {print join("\t", @F);}' ${TARGETDIR}/tmp/temp${NUM}.edge_broden.bed | sort -k 4 - > ${TARGETDIR}/tmp/temp${NUM}.edge_broden.sorted.bed
check_error $?

echo "perl summarizeEdge.pl ${TARGETDIR}/tmp/temp${NUM}.edge.sorted.bed > ${TARGETDIR}/tmp/temp${NUM}.edge2count.bed"
perl summarizeEdge.pl ${TARGETDIR}/tmp/temp${NUM}.edge.sorted.bed > ${TARGETDIR}/tmp/temp${NUM}.edge2count.bed
check_error $?

echo "perl summarizeEdge_broden.pl ${TARGETDIR}/tmp/temp${NUM}.edge_broden.sorted.bed > ${TARGETDIR}/tmp/temp${NUM}.edge_broden2count.bed"
perl summarizeEdge_broden.pl ${TARGETDIR}/tmp/temp${NUM}.edge_broden.sorted.bed > ${TARGETDIR}/tmp/temp${NUM}.edge_broden2count.bed
check_error $?

echo "perl mergeEdges.pl ${TARGETDIR}/tmp/temp${NUM}.edge2count.bed ${TARGETDIR}/tmp/temp${NUM}.edge_broden2count.bed > ${TARGETDIR}/tmp/temp${NUM}.edge2count_merge.bed"
perl mergeEdges.pl ${TARGETDIR}/tmp/temp${NUM}.edge2count.bed ${TARGETDIR}/tmp/temp${NUM}.edge_broden2count.bed > ${TARGETDIR}/tmp/temp${NUM}.edge2count_merge.bed
check_error $?

rm -rf ${TARGETDIR}/tmp/temp${NUM}.bed12
rm -rf ${TARGETDIR}/tmp/temp${NUM}.edge.bed
rm -rf ${TARGETDIR}/tmp/temp${NUM}.edge.sorted.bed
rm -rf ${TARGETDIR}/tmp/temp${NUM}.edge_broden.bed
rm -rf ${TARGETDIR}/tmp/temp${NUM}.edge_broden.sorted.bed
rm -rf ${TARGETDIR}/tmp/temp${NUM}.edge2count.bed
rm -rf ${TARGETDIR}/tmp/temp${NUM}.edge_broden2count.bed

##########

echo "R --vanilla --slave --args ${TARGETDIR}/tmp/temp${NUM}.edge2count_merge.bed ${TARGETDIR}/tmp/temp${NUM}.edge2count_bayes.bed < bayesBetaFilter.R"
R --vanilla --slave --args ${TARGETDIR}/tmp/temp${NUM}.edge2count_merge.bed ${TARGETDIR}/tmp/temp${NUM}.edge2count_bayes.bed < bayesBetaFilter.R

# _COMMENT_OUT_


