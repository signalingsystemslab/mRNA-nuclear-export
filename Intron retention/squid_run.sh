#!/bin/bash

# source venv/bin/activate

echo $date

ANNOT="/opt/ngs_indexes/models/mm/mm10/gencode.vM14.annotation.gtf"
IN_DIR="/home/diane/merged_bams_mm10_vM14"

N_THREADS="8"

export ANNOT OUT_DIR IN_DIR N_THREADS
 
run_squid(){
    BAM_DIR=$1
	echo $BAM_DIR
    BAMS=$(ls $IN_DIR*/$BAM_DIR/*[0,5].rep[1-3].bam | sed -e "s@^\(${IN_DIR}.*\)/\(.*\)\$@\2#\1@"| sort | sed 's@^\(.*\)#\(.*\)$@\2/\1@' | tr '\n' ',' | sed 's@,$@@') 
    echo $BAMS	
    # OUT_NAME=$(echo $BAM_DIR | cut -d"/" -f 5)
    
    # echo $OUT_NAME
    echo $ANNOT
    echo $OUT_DIR
    

    mkdir -p "$OUT_DIR/$BAM_DIR"

    python SQUID.py --GTF $ANNOT --align $BAMS --output $OUT_DIR/$BAM_DIR --lib first --read P --Cal count --p $N_THREADS --resume true

}

export -f run_squid

OUT_DIR="SQUID/Naive"
mkdir -p $OUT_DIR
parallel -j $N_THREADS --no-notice run_squid {} ::: $(ls $IN_DIR/ | grep Nuc.Naive.rep[1-3])
gzip -9 $OUT_DIR/*/Result/intron_PI.txt

OUT_DIR="SQUID/LPA"
mkdir -p $OUT_DIR
parallel -j $N_THREADS --no-notice run_squid {} ::: $(ls $IN_DIR/ | grep Nuc.LPA.rep[1-3])
gzip -9 $OUT_DIR/*/Result/intron_PI.txt

unset -f run_squid

date

unset ANNOT OUT_DIR IN_DIR N_THREADS

# deactivate

