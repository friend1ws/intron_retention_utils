#! /bin/sh
#$ -S /bin/sh
#$ -cwd


source ../lib/config.sh
source ../lib/utility.sh
PATH=${BLAT_PATH}:$PATH
referece=${BLAT_REF}

PATH=${BAEDTOOLS_PATH}:${PATH}
PATH=${SAMTOOLS_PATH}:${PATH}

SAMPLE=$1

SEQDIR=${OUTPUTDIR}/${SAMPLE}/bam
TARGETDIR=${OUTPUTDIR}/${SAMPLE}/intron


if [ -f ${SEQDIR}/sequence.bam ]
then
    echo "${SEQDIR}/sequence.bam exists!"
else
    echo "${SEQDIR}/sequence.bam does not exist!"
    exit
fi

CURLOGDIR=${LOGDIR}/${SAMPLE}
LOGSTR=-e\ ${CURLOGDIR}\ -o\ ${CURLOGDIR}

echo ${LOGSTR}

if [ -d ${TARGETDIR}/tmp ]
then
    echo "${TARGETDIR}/tmp exists."
else
    echo "${TARGETDIR}/tmp does not exits."
    mkdir -p ${TARGETDIR}/tmp
fi

if [ -d ${LOGDIR} ]
then
    echo "${LOGDIR} exists."
else
    echo "${LOGDIR} does not exits."
    mkdir -p ${LOGDIR}
fi


FILECOUNT=`find ${DBDIR}/interval_list_hg19_nongap/*.interval_list | wc -l`

job_getJunction=getIntron.${SAMPLE}

echo "qsub -t 1-${FILECOUNT}:1 -l s_vmem=8G,mem_req=8 -N ${job_getJunction} ${LOGSTR} getIntron_part.sh ${SAMPLE}"
qsub -t 1-${FILECOUNT}:1 -l s_vmem=8G,mem_req=8 -N ${job_getJunction} ${LOGSTR} getIntron_part.sh ${SAMPLE}



job_catIntron=catIntron.${SAMPLE}
echo "qsub -l s_vmem=8G,mem_req=8 -N ${job_catIntron} -hold_jid ${job_getJunction} ${LOGSTR} catIntron.sh ${SAMPLE}"
qsub -l s_vmem=8G,mem_req=8 -N ${job_catIntron} -hold_jid ${job_getJunction} ${LOGSTR} catIntron.sh ${SAMPLE}



