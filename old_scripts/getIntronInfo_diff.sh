#! /bin/sh
#$ -S /bin/sh
#$ -cwd


source ../lib/config.sh
source ../lib/utility.sh
PATH=${BLAT_PATH}:$PATH
referece=${BLAT_REF}

PATH=${BAEDTOOLS_PATH}:${PATH}
PATH=${SAMTOOLS_PATH}:${PATH}

SAMPLE1=$1
SAMPLE2=$2
# THRES=$2

TARGETDIR=${RESULTDIR}/${SAMPLE1}-${SAMPLE2}/intron


CURLOGDIR=${LOGDIR}/${SAMPLE1}-${SAMPLE2}
LOGSTR=-e\ ${CURLOGDIR}\ -o\ ${CURLOGDIR}

echo ${LOGSTR}

if [ -d ${TARGETDIR}/tmp ]
then
    echo "${TARGETDIR}/tmp exists."
else
    echo "${TARGETDIR}/tmp does not exits."
    mkdir -p ${TARGETDIR}/tmp
fi


FILECOUNT=`find ${DBDIR}/interval_list_hg19_nongap/*.interval_list | wc -l`

job_getJunction_diff=getJunction_diff.${SAMPLE1}-${SAMPLE2}

echo "qsub -t 1-${FILECOUNT}:1 -l s_vmem=8G,mem_req=8 -N ${job_getJunction_diff} ${LOGSTR} getIntron_diff.sh ${SAMPLE1} ${SAMPLE2}"
qsub -t 1-${FILECOUNT}:1 -l s_vmem=8G,mem_req=8 -N ${job_getJunction_diff} ${LOGSTR} getIntron_diff.sh ${SAMPLE1} ${SAMPLE2}



job_catCandSplicing_diff=catCandSplicing_diff.${SAMPLE1}-${SAMPLE2}2
echo "qsub -l s_vmem=8G,mem_req=8 -N ${job_catCandSplicing_diff} -hold_jid ${job_getJunction_diff} ${LOGSTR} catCandSplicing_diff.sh ${SAMPLE1} ${SAMPLE2}"
# qsub -l s_vmem=8G,mem_req=8 -N ${job_catCandSplicing_diff} -hold_jid ${job_getJunction_diff} ${LOGSTR} catCandSplicing_diff.sh ${SAMPLE1} ${SAMPLE2}



