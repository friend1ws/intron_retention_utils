#! /bin/sh
#$ -S /bin/sh
#$ -cwd

source ../lib/config.sh
source ../lib/utility.sh
PATH=${BLAT_PATH}:$PATH
referece=${BLAT_REF}

PATH=${BEDTOOLS_PATH}:${PATH}
PATH=${SAMTOOLS_PATH}:${PATH}

SAMPLE1=$1
SAMPLE2=$2
INTRON_REF=$3
# NUM=${SGE_TASK_ID}
# NUM=3
echo $R_LIBS
export R_LIBS=/home/yshira/.R
echo $R_LIBS

INPUTDIR1=${OUTPUTDIR}/${SAMPLE1}/intron
INPUTDIR2=${OUTPUTDIR}/${SAMPLE2}/intron
TARGETDIR=${RESULTDIR}/${SAMPLE1}-${SAMPLE2}/intron

MAPG_REF=/home/yshira/RNAseq/script/intron_ref/mappedBaseCount_merged.txt
MAPG_T=/home/yshira/RNAseq/data/output/${SAMPLE1}/expression/mappedBaseCount.txt
# INTERVAL=${DBDIR}/interval_list_hg19_nongap

# sleep 
# sh sleep.sh
 
check_mkdir ${TARGETDIR}
##########
# :<<'__COMMENT_OUT__'

echo "perl mergeTwoIntron.pl ${INPUTDIR1}/intron2count.bed ${INPUTDIR2}/intron2count.bed > ${TARGETDIR}/intron2count.bed"
perl mergeTwoIntron.pl ${INPUTDIR1}/intron2count.bed ${INPUTDIR2}/intron2count.bed > ${TARGETDIR}/intron2count.bed
check_error $?

# __COMMENT_OUT__
echo "perl addRef.pl ${TARGETDIR}/intron2count.bed ${INTRON_REF} > ${TARGETDIR}/intron2count_ref.bed"
perl addRef.pl ${TARGETDIR}/intron2count.bed ${INTRON_REF} > ${TARGETDIR}/intron2count_ref.bed
check_error $?

# __COMMENT_OUT__

echo "R --vanilla --slave --args ${TARGETDIR}/intron2count_ref.bed ${MAPG_REF} ${MAPG_T} ${TARGETDIR}/intron2count_eb.bed < proc_EBcall.R"
R --vanilla --slave --args ${TARGETDIR}/intron2count_ref.bed ${MAPG_REF} ${MAPG_T} ${TARGETDIR}/intron2count_eb.bed < proc_EBcall.R
check_error $?

##########
:<<_COMMENT_OUT_
# obtain the coverage arround the junction start sites

echo "intersectBed -a ${INPUTDIR1}/tmp/temp${NUM}.bed12 -b ${TARGETDIR}/tmp/temp${NUM}.edge2count_merge.bed -split -wao > ${TARGETDIR}/tmp/temp${NUM}.edge_1.bed"
intersectBed -a ${INPUTDIR1}/tmp/temp${NUM}.bed12 -b ${TARGETDIR}/tmp/temp${NUM}.edge2count_merge.tmp.bed -split -wao > ${TARGETDIR}/tmp/temp${NUM}.edge_1.bed
check_error $?

echo "intersectBed -a ${INPUTDIR2}/tmp/temp${NUM}.bed12 -b ${TARGETDIR}/tmp/temp${NUM}.edge2count_merge.bed -split -wao > ${TARGETDIR}/tmp/temp${NUM}.edge_2.bed"
intersectBed -a ${INPUTDIR2}/tmp/temp${NUM}.bed12 -b ${TARGETDIR}/tmp/temp${NUM}.edge2count_merge.tmp.bed -split -wao > ${TARGETDIR}/tmp/temp${NUM}.edge_2.bed
check_error $?

echo "perl summarizeEdge_diff.pl ${TARGETDIR}/tmp/temp${NUM}.edge_1.bed > ${TARGETDIR}/tmp/temp${NUM}.edge2count_1.bed"
perl summarizeEdge_diff.pl ${TARGETDIR}/tmp/temp${NUM}.edge_1.bed > ${TARGETDIR}/tmp/temp${NUM}.edge2count_1.bed
check_error $?

echo "perl summarizeEdge_diff.pl ${TARGETDIR}/tmp/temp${NUM}.edge_2.bed > ${TARGETDIR}/tmp/temp${NUM}.edge2count_2.bed"
perl summarizeEdge_diff.pl ${TARGETDIR}/tmp/temp${NUM}.edge_2.bed > ${TARGETDIR}/tmp/temp${NUM}.edge2count_2.bed
check_error $?

echo "perl summarizeTwoIntron.pl ${TARGETDIR}/tmp/temp${NUM}.edge2count_1.bed ${TARGETDIR}/tmp/temp${NUM}.edge2count_2.bed ${TARGETDIR}/tmp/temp${NUM}.edge2count_merge.tmp.bed > ${TARGETDIR}/tmp/temp${NUM}.edge2count_merge.bed"
perl summarizeTwoIntron.pl ${TARGETDIR}/tmp/temp${NUM}.edge2count_1.bed ${TARGETDIR}/tmp/temp${NUM}.edge2count_2.bed ${TARGETDIR}/tmp/temp${NUM}.edge2count_merge.tmp.bed > ${TARGETDIR}/tmp/temp${NUM}.edge2count_merge.bed
check_error $?

##########

echo "R --vanilla --slave --args ${TARGETDIR}/tmp/temp${NUM}.edge2count_merge.bed ${TARGETDIR}/tmp/temp${NUM}.edge2count_fis.bed < proc_fisherTest.R"
R --vanilla --slave --args ${TARGETDIR}/tmp/temp${NUM}.edge2count_merge.bed ${TARGETDIR}/tmp/temp${NUM}.edge2count_fis.bed < proc_fisherTest.R
check_error $?

_COMMENT_OUT_

