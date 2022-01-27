#!/bin/bash

#  03_trimming.sh
#  
#
#  Created by Diane on 10/26/16.
#

N_THREADS="24"
# take input from following directory
export IN_DIR=$1
# write output files to following directory
export OUT_DIR="03_trimmed_index_universal/"
# create OUT_DIR if it doesn't exist
mkdir -p $OUT_DIR

export SUFF=".fq"
LOG="03_trimming_index_universal"

mkdir -p ${OUT_DIR}cutadapt

trim_function () { 
    FASTQ_FILE=$1
    BASE=$(basename $FASTQ_FILE | sed 's@\(\.\|_\)1\(\.fq\|\.fastq\)@@' | cut -d"/" -f2) 
    MATE=$(echo $FASTQ_FILE | sed 's@1\.f@2\.f@')

    date > $OUT_DIR$BASE.tmp.log
    pwd >> $OUT_DIR$BASE.tmp.log
    
    echo "-----------------" >> $OUT_DIR$BASE.tmp.log 
    echo "Trimming adapter for $BASE" >> $OUT_DIR$BASE.tmp.log
#     echo "cutadapt -f fastq -e 0.1 -O 8 -a index=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A universal=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -m 23 -o ${OUT_DIR}${BASE}_1${SUFF} -p ${OUT_DIR}${BASE}_2${SUFF} ${IN_DIR}${BASE}_1${SUFF} ${IN_DIR}${BASE}_2${SUFF}" >> $OUT_DIR$BASE.tmp.log
#    echo "cutadapt -f fastq -e 0.1 -O 8 -a common=AGATCGGAAGAGC -A common=AGATCGGAAGAGC -m 23 -o ${OUT_DIR}${BASE}_1${SUFF} -p ${OUT_DIR}${BASE}_2${SUFF} ${IN_DIR}${BASE}_1${SUFF} ${IN_DIR}${BASE}_2${SUFF}" >> $OUT_DIR$BASE.tmp.log
    (cutadapt -f fastq \
        -e 0.1 \
        -O 8 \
        -a index=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        -A universal=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
        -m 23 \
        -o ${OUT_DIR}${BASE}_1${SUFF} \
        -p ${OUT_DIR}${BASE}_2${SUFF} \
        $FASTQ_FILE \
        $MATE 2>&1 & echo $! > pid) | tee -a ${OUT_DIR}cutadapt/$BASE.tmp.log &

#        -a index=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
#        -A universal=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
#	    -a common=AGATCGGAAGAGC \
#        -A common=AGATCGGAAGAGC \
#    if [$(cat pid | wc -l) -ne 0]; then
#        echo $(cat pid)
#        memory_script.sh $(cat pid) 5 &> $OUT_DIR$BASE.mem 
#    fi

    wait
}


# for FASTQ_FILE in $(find ${IN_DIR} -name "*1${SUFF}" -maxdepth 1)
# do
# 	trims $FASTQ_FILE
# done


export -f trim_function 
parallel -j $N_THREADS --no-notice trim_function ::: $(ls ${IN_DIR}*1.f*q) 
unset -f trim_function

cat ${OUT_DIR}*.tmp.log > $LOG.log
rm ${OUT_DIR}*.tmp.log

#cat ${OUT_DIR}*.mem > $LOG.mem
#rm ${OUT_DIR}*mem
 
echo -e "\n----------------" >> $LOG.log
echo "Creating FastQC reports" >> $LOG.log

mkdir -p ${OUT_DIR}fastqc
echo "fastqc -t $N_THREADS --outdir ${OUT_DIR}fastqc ${OUT_DIR}*${SUFF}" >> $LOG.log
(fastqc -t $N_THREADS --outdir ${OUT_DIR}fastqc ${OUT_DIR}*${SUFF} 2>&1 & echo $! >pid)| tee -a $LOG.log &
# memory_script.sh $(cat pid) 1 &>> $LOG.mem
wait

echo "----------" >> $LOG.log
echo "Creating MultiQC reports" >> $LOG.log
mkdir -p ${OUT_DIR}multiqc
echo "multiqc $OUT_DIR  --outdir ${OUT_DIR}multiqc --interacive" >> $LOG.log
multiqc $OUT_DIR --outdir ${OUT_DIR}multiqc --interactive -f 2>&1 | tee -a $LOG.log 

date >> $LOG.log

echo "Compressing all fastq files from $IN_DIR" >> $LOG.log
# echo "parallel -j $N_THREADS gzip -9 ::: ${IN_DIR}*${SUFF}" >> $LOG.log
parallel -j $N_THREADS --no-notice gzip -9 ::: $( ls ${IN_DIR}*.f*q) 2>&1 | tee -a $LOG.log

date >> $LOG.log

unset IN_DIR OUT_DIR SUFF
rm pid
