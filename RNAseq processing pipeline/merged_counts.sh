#!/bin/sh

#  08_counts.sh
#  
#
#  Created by Diane on 10/31/16.
#

N_THREADS=24

# take input from following directory
IN_DIR="merged_bams_mm10_vM14"
# write output files to following directory
OUT_DIR="merged_counts"
mkdir -p $OUT_DIR

# Genome annotation
GENOME_ANNOT="/opt/ngs_indexes/models/mm/mm10/gencode.vM14.annotation.gtf"

for DIR2 in $(ls -d ${IN_DIR}*/*/ | xargs -d '\n' -n 1 basename | sort | uniq )
do 
    echo $DIR2
    mkdir -p ${OUT_DIR}/${DIR2}"/D40000"
    mkdir -p ${OUT_DIR}/${DIR2}"/D40000/exons"
    mkdir -p ${OUT_DIR}/${DIR2}"/D40000/all"

    featureCounts -T 8 -s 2 -S fr -p -B -C -P -d 40 -D 40000 -t exon -g gene_id -a ${GENOME_ANNOT} -o ${OUT_DIR}/${DIR2}/D40000/exons/exon_counts.txt ${IN_DIR}*/${DIR2}/*.bam 2>&1 &
    
    wait
    
    featureCounts -T 8 -s 2 -S fr -p -B -C -P -d 40 -D 40000 -t gene -g gene_id -a ${GENOME_ANNOT} -o ${OUT_DIR}/${DIR2}/D40000/all/all_counts.txt ${IN_DIR}*/${DIR2}/*.bam 2>&1 &
    
    wait
    
done
