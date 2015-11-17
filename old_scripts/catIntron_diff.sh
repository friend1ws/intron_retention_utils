#! /bin/sh
#$ -S /bin/sh
#$ -cwd

source ../lib/config.sh
source ../lib/utility.sh

PATH=${BLAT_PATH}:$PATH
referece=${BLAT_REF}


SAMPLE1=$1
SAMPLE2=$2

TARGETDIR=${RESULTDIR}/${SAMPLE1}-${SAMPLE2}/intron

PATH=${BEDTOOLS_PATH}:${PATH}
PATH=${SAMTOOLS_PATH}:${PATH}



FILECOUNT=`find ${DBDIR}/interval_list_hg19_nongap/*.interval_list | wc -l`

echo -n > ${TARGETDIR}/intron2info.bed
for i in `seq 1 1 ${FILECOUNT}`
do
    echo "cat ${TARGETDIR}/tmp/temp${i}.edge2count_fis.bed >> ${TARGETDIR}/intron2info.bed"
    cat ${TARGETDIR}/tmp/temp${i}.edge2count_fis.bed >> ${TARGETDIR}/intron2info.bed
    check_error $?
done

