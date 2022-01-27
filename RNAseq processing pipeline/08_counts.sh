#!/bin/sh

#  08_counts.sh
#  
#
#  Created by Diane on 10/31/16.
#

N_THREADS=24

# take input from following directory
IN_DIR="06_filtering/b_deRibo/"
# write output files to following directory
OUT_DIR="08_counts/"

LOG="08_counts"

# create OUT_DIR if it doesn't exist
mkdir -p $OUT_DIR
mkdir -p "${OUT_DIR}all"
mkdir -p "${OUT_DIR}exons"

# Genome annotation
GENOME_ANNOT="/opt/ngs_indexes/models/mm/mm10/gencode.vM14.annotation.gtf"

date >> $LOG.log

(featureCounts -T $N_THREADS -s 2 -S fr -p -B -C -P -d 40 -D 4000 -t gene -g gene_id -a ${GENOME_ANNOT} -o ${OUT_DIR}all/all_counts.txt $IN_DIR*[^{.junk},^{.ribo}].bam 2>&1 & echo $! > pid) | tee -a $LOG.log
# memory_script.sh $(cat pid) 1 >> $LOG.mem

wait

cat $OUT_DIR/all/all_counts.txt.summary | sed 's@\.bam@_all\.bam@g' > $OUT_DIR/all/all_counts_2.txt.summary
mv $OUT_DIR/all/all_counts_2.txt.summary $OUT_DIR/all/all_counts.txt.summary

(featureCounts -T $N_THREADS -s 2 -S fr -p -B -C -P -d 40 -D 4000 -t exon -g gene_id -a ${GENOME_ANNOT} -o ${OUT_DIR}exons/exon_counts.txt $IN_DIR*[^{.junk},^{.ribo}].bam 2>&1 & echo $! > pid) | tee -a $LOG.log
# memory_script.sh $(cat pid) 1 >> $LOG.mem

wait

cat $OUT_DIR/exons/exon_counts.txt.summary | sed 's@\.bam@_exons\.bam@g' > $OUT_DIR/exons/exon_counts_2.txt.summary
mv $OUT_DIR/exons/exon_counts_2.txt.summary $OUT_DIR/exons/exon_counts.txt.summary

mkdir -p ${OUT_DIR}multiqc
echo "multiqc ${OUT_DIR}all/  -f --outdir ${OUT_DIR}multiqc_all --interacive" >> $LOG.log
multiqc ${OUT_DIR} -f --outdir ${OUT_DIR}multiqc --interactive 2>&1 | tee -a $LOG.log

wait

# mkdir -p ${OUT_DIR}multiqc_exons
# echo "multiqc ${OUT_DIR}exons/  -f --outdir ${OUT_DIR}multiqc_exons --interacive" >> $LOG.log         date >> $LOG.log
# multiqc ${OUT_DIR}exons/ -f --outdir ${OUT_DIR}multiqc_exons --interactive 2>&1 | tee -a $LOG.log
# wait

rm pid
