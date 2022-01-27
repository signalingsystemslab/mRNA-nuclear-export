#!/bin/sh

#  02_consolidated_fastq_michael.sh
#  
#
#  Created by Diane on 10/25/16.
#

date

if ls 00_qseq/*.gz 2>/dev/null
then
    parallel -j 24 --no-notice gzip -d ::: 00_qseq/*.gz
fi
wait

NOW=$(date +"%Y_%m_%d")
OUT_LOG="02_consolidated_$NOW.log"
# OUT_MEM="02_consolidated_$NOW.mem"

date > $OUT_LOG 
pwd >> $OUT_LOG

# take input from following directory, with following suffix
IN_DIR="00_qseq"
# write output files to following directory, with following suffixp
OUT_DIR="02_consolidated/"

N_THREADS="8"
BARCODE="barcodes.txt"

# create OUT_DIR if it doesn't exist
mkdir -p $OUT_DIR

# convert all qseq files to fastq and demultiplex
laneRunner.py -i $IN_DIR -o $OUT_DIR -b $BARCODE -p $N_THREADS -s 2>&1 >> $OUT_LOG &
# memory_script.sh $! 10 > $OUT_MEM
wait

date >> $OUT_LOG

# ./run_fastqscreen.sh

# echo "/n----------" | tee -a $OUT_LOG $OUT_MEM
# echo "Creating FastQC reports" | tee -a $OUT_LOG $OUT_MEM
# mkdir -p ${OUT_DIR}fastqc
# echo "fastqc -t $N_THREADS --outdir ${OUT_DIR}fastqc ${OUT_DIR}*.fastq" >> $OUT_LOG
# fastqc -t $N_THREADS --outdir ${OUT_DIR}fastqc ${OUT_DIR}*.fastq 2>&1 >> $OUT_LOG &
# memory_script.sh $! 10 >> $OUT_MEM

echo "Compressing all qseq files from $IN_DIR" >> $OUT_LOG
parallel -j $N_THREADS --no-notice gzip -9 {} ::: ${IN_DIR}/*qseq.txt 2>&1 >> $OUT_LOG 

rm -r ${OUT_DIR}parts
rm -r ${OUT_DIR}.concat*
rm -r ${OUT_DIR}.demux*

# remove lane number in the name
for f in ${OUT_DIR}*.f*q
do 
    mv $f $(echo $f | sed 's@\..\.@\.@')
done


