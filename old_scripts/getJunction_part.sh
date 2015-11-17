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
TARGETDIR=${OUTPUTDIR}/${SAMPLE}/junction

INTERVAL=${DBDIR}/interval_list_hg19_nongap

# sleep 
# sh sleep.sh
 

REGION_A=`head -n1 ${INTERVAL}/${NUM}.interval_list | awk '{split($0, ARRAY, "-"); print ARRAY[1]}'`
REGION_B=`tail -n1 ${INTERVAL}/${NUM}.interval_list | awk '{split($0, ARRAY, "-"); print ARRAY[2]}'`
REGION="${REGION_A}-${REGION_B}"
echo ${REGION} 

##########
# : <<'_COMMENT_OUT__'

echo "samtools view -h ${SEQDIR}/sequence.bam ${REGION} >  ${TARGETDIR}/tmp/temp${NUM}.sam"
samtools view -h ${SEQDIR}/sequence.bam ${REGION} > ${TARGETDIR}/tmp/temp${NUM}.sam
check_error $?


echo "perl filterRead.pl ${TARGETDIR}/tmp/temp${NUM}.sam 60 > ${TARGETDIR}/tmp/temp${NUM}.sam.filt"
perl filterRead.pl ${TARGETDIR}/tmp/temp${NUM}.sam 60 > ${TARGETDIR}/tmp/temp${NUM}.sam.filt
check_error $?

echo "samtools view -h -bS ${TARGETDIR}/tmp/temp${NUM}.sam.filt > ${TARGETDIR}/tmp/temp${NUM}.bam"
samtools view -h -bS ${TARGETDIR}/tmp/temp${NUM}.sam.filt > ${TARGETDIR}/tmp/temp${NUM}.bam
check_error $?


echo "bamToBed -i ${TARGETDIR}/tmp/temp${NUM}.bam -bed12 > ${TARGETDIR}/tmp/temp${NUM}.bed12"
bamToBed -i ${TARGETDIR}/tmp/temp${NUM}.bam -bed12 > ${TARGETDIR}/tmp/temp${NUM}.bed12
check_error $?


echo "perl getSpliceJunction.pl ${TARGETDIR}/tmp/temp${NUM}.bed12 ${THRES}  > ${TARGETDIR}/tmp/junction2count${NUM}.bed"
perl getSpliceJunction.pl ${TARGETDIR}/tmp/temp${NUM}.bed12 ${THRES}  > ${TARGETDIR}/tmp/junction2count${NUM}.bed
check_error $?


##########
# obtain the coverage arround the junction start sites
perl -F"\t" -a -n -e 'if ($F[3] >= 3) {print $F[0] . "\t" . ($F[1] - 1) . "\t" . $F[1] . "\t" . join("\t", @F);}' ${TARGETDIR}/tmp/junction2count${NUM}.bed > ${TARGETDIR}/tmp/junction2start${NUM}.bed
check_error $?

perl -F"\t" -a -n -e 'if ($F[3] >= 3) {print $F[0] . "\t" . $F[2] . "\t" . ($F[2] + 1) . "\t" . join("\t", @F);}' ${TARGETDIR}/tmp/junction2count${NUM}.bed > ${TARGETDIR}/tmp/junction2end${NUM}.bed
check_error $?

# echo "coverageBed -a ${TARGETDIR}/tmp/temp${NUM}.bed12 -b ${TARGETDIR}/tmp/junction2start${NUM}.bed -split -d > ${TARGETDIR}/tmp/junction2countCoverage${NUM}.tmp.bed"
# coverageBed -a ${TARGETDIR}/tmp/temp${NUM}.bed12 -b ${TARGETDIR}/tmp/junction2start${NUM}.bed -split -d > ${TARGETDIR}/tmp/junction2countCoverage${NUM}.tmp.bed
# check_error $?

echo "coverageBed -a ${TARGETDIR}/tmp/temp${NUM}.bam -b ${TARGETDIR}/tmp/junction2start${NUM}.bed -split -d > ${TARGETDIR}/tmp/junction2countCoverage${NUM}.start.tmp.bed"
coverageBed -abam ${TARGETDIR}/tmp/temp${NUM}.bam -b ${TARGETDIR}/tmp/junction2start${NUM}.bed -split -d > ${TARGETDIR}/tmp/junction2countCoverage${NUM}.start.tmp.bed
check_error $?

echo "coverageBed -a ${TARGETDIR}/tmp/temp${NUM}.bam -b ${TARGETDIR}/tmp/junction2end${NUM}.bed -split -d > ${TARGETDIR}/tmp/junction2countCoverage${NUM}.end.tmp.bed"
coverageBed -abam ${TARGETDIR}/tmp/temp${NUM}.bam -b ${TARGETDIR}/tmp/junction2end${NUM}.bed -split -d > ${TARGETDIR}/tmp/junction2countCoverage${NUM}.end.tmp.bed
check_error $?


# _COMMENT_OUT__

# perl -F"\t" -a -n -e 'print join("\t", @F[(3, 4, 5, 6, 8)]);' ${TARGETDIR}/tmp/junction2countCoverage${NUM}.tmp.bed > ${TARGETDIR}/tmp/junction2countCoverage${NUM}.tmp2.bed
# check_error $?

echo "perl mergeStartEndCoverage.pl ${TARGETDIR}/tmp/junction2countCoverage${NUM}.start.tmp.bed ${TARGETDIR}/tmp/junction2countCoverage${NUM}.end.tmp.bed > ${TARGETDIR}/tmp/junction2countCoverage${NUM}.tmp2.bed"
perl mergeStartEndCoverage.pl ${TARGETDIR}/tmp/junction2countCoverage${NUM}.start.tmp.bed ${TARGETDIR}/tmp/junction2countCoverage${NUM}.end.tmp.bed > ${TARGETDIR}/tmp/junction2countCoverage${NUM}.tmp2.bed

echo "R --vanilla --slave --args ${TARGETDIR}/tmp/junction2countCoverage${NUM}.tmp2.bed ${TARGETDIR}/tmp/junction2countCoverage${NUM}.bed < bayesBetaFilter.R"
R --vanilla --slave --args ${TARGETDIR}/tmp/junction2countCoverage${NUM}.tmp2.bed ${TARGETDIR}/tmp/junction2countCoverage${NUM}.bed < bayesBetaFilter.R
check_error $?

##########




