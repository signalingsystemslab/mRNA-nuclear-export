#!/bin/bash

N_THREADS=24

export CHROM="/opt/ngs_indexes/star/mm10_no_alt_analysis_set_ENCODE.gencode_vM14_annotation.101bp/chrNameLength.txt"
export IN_DIR="06_filtering/b_deRibo/"
export OUT_DIR="07_tracks_no_multi_hit/"

LOG="07_tracks_unique"

date > $LOG.log

mkdir -p $OUT_DIR

track_function () {
BAM_FILE=$1
NAME=$(echo $BAM_FILE | cut -d"/" -f 3 | cut -d"." -f 1)

echo "" >> $NAME.tmp.log 
echo "-----------------" >> $NAME.tmp.log 
date >> $NAME.tmp.log
echo "Making tracks for $NAME" >> $NAME.tmp.log
(bam2wig.py -i "$BAM_FILE" -s $CHROM -o ${OUT_DIR}${NAME} -d '1+-,1-+,2++,2--' -t 1000000000 -u 2>&1 & echo $! >pid) | tee -a $NAME.tmp.log &
# memory_script.sh $(cat pid) 1 >> $NAME.tmp.mem
wait 

# (wigToBigWig $OUT_DIR${NAME}.Forward.wig $CHROM $OUT_DIR${NAME}.Forward.bw 2>&1 & echo $! >pid) | tee -a $NAME.tmp.log &
# # memory_script.sh $(cat pid) 1 >> $NAME.tmp.mem
# wait
 
# (wigToBigWig $OUT_DIR${NAME}.Reverse.wig $CHROM $OUT_DIR${NAME}.Reverse.bw 2>&1 & echo $! >pid) | tee -a $NAME.tmp.log &
# # memory_script.sh $(cat pid) 1 >> $NAME.tmp.mem
# wait

echo "" >> $NAME.tmp.log
date >> $NAME.tmp.log

}

export -f track_function
parallel -j $N_THREADS --no-notice track_function ::: $(ls $IN_DIR*[^{.junk},^{.ribo}].bam)
unset -f track_function

rm $OUT_DIR*.wig

cat *.tmp.log >> $LOG.log
rm *.tmp.log

# cat *.tmp.mem >> $LOG.mem
# rm *.tmp.mem
 
unset CHROM IN_DIR OUT_DIR
