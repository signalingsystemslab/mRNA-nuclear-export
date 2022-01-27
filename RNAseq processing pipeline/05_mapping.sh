#!/bin/bash

#  05_alignment.sh
#  
#
#  Created by Diane on 10/26/16.
#

N_THREADS="12"

# take input from following directory
export IN_DIR="03_trimmed_index_universal/"
# write output files to following directory
export OUT_DIR="05_alignment/"
# create OUT_DIR
mkdir -p $OUT_DIR

LOG="05_alignment"

# ref genome
export GENOME_DIR="/opt/ngs_indexes/star/mm10_no_alt_analysis_set_ENCODE.gencode_vM14_annotation.101bp"
export GENOME_ANNOT_BED="/opt/ngs_indexes/models/mm/mm10/gencode.vM14.annotation.bed"
export GENOME_ANNOT_GTF="/opt/ngs_indexes/models/mm/mm10/gencode.vM14.annotation.gtf"
# export GENOME_DIR_RSEM="/opt/ngs_indexes/rsem/mm10"

# file suffix
export SUFF=".fq"

date > $LOG.log 
pwd >> $LOG.log

STAR --genomeDir $GENOME_DIR --genomeLoad LoadAndExit

wait 

map_function () {
    FASTQ_FILE=$1
    BASE=$(basename -s "_1$SUFF" $FASTQ_FILE)
    mkdir -p ${OUT_DIR}${BASE}
    
    date > $OUT_DIR$BASE.tmp.log
    pwd >> $OUT_DIR$BASE.tmp.log
    
    echo "----------------" >> $OUT_DIR$BASE.tmp.log 
    echo "Mapping $BASE to ref genome $GENOME_DIR" >> $OUT_DIR$BASE.tmp.log
    (STAR --genomeDir $GENOME_DIR \
         --runThreadN 2 \
         --genomeLoad LoadAndKeep \
         --readFilesIn ${IN_DIR}${BASE}_1${SUFF} ${IN_DIR}${BASE}_2${SUFF} \
         --outSAMunmapped Within \
         --outSAMtype BAM SortedByCoordinate \
         --limitBAMsortRAM 17179869184 \
         --outFilterType BySJout \
         --outFilterMultimapNmax 50 \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --alignIntronMin 20 \
         --alignIntronMax 1000000 \
         --alignMatesGapMax 1000000 \
         --outFilterMismatchNoverLmax 0.04 \
         --seedSearchStartLmax 30 \
         --alignEndsType EndToEnd \
         --outFileNamePrefix "$OUT_DIR${BASE}/" 2>&1 & echo $! >pid)| tee -a $OUT_DIR$BASE.tmp.log &
    
#    if cat pid 2>/dev/null 
#    then
#        memory_script.sh $(cat pid) 5 &> $OUT_DIR$BASE.tmp.mem
#    fi

    wait
    mv "$OUT_DIR${BASE}/Aligned.sortedByCoord.out.bam" "${OUT_DIR}${BASE}/$BASE.bam"
    # mv "$OUT_DIR${BASE}/Aligned.toTranscriptome.out.bam" "${OUT_DIR}${BASE}/$BASE.toTranscriptome.bam"
}
 
export -f map_function
parallel -j $N_THREADS --no-notice map_function ::: ${IN_DIR}*1${SUFF}
unset -f map_function
   
STAR --genomeDir $GENOME_DIR --genomeLoad Remove
 
cat ${OUT_DIR}*.tmp.log > $LOG.log
rm ${OUT_DIR}*.tmp.log

#echo "-------------" > $LOG.mem
#echo "STAR" >> $LOG.mem
#cat ${OUT_DIR}*.tmp.mem >> $LOG.mem
#rm ${OUT_DIR}*.tmp.mem

wait

echo "--------------------" >> $LOG.log
echo "Index BAM files" >> $LOG.log
# echo "parallel -j $N_THREADS --no-notice samtools index ::: $OUT_DIR*/*.bam" >> $LOG.log
parallel -j $(echo `expr $N_THREADS \* 2`) --no-notice '(echo {}; samtools index {})' ::: $(ls $OUT_DIR*/*[^e].bam) 2>&1 | tee -a $LOG.log
 
echo "--------------------" >> $LOG.log
echo "Infer details of experiments" >> $LOG.log
# echo "parallel -j $N_THREADS --no-notice infer_experiment.py -i {} -r $GENOME_ANNOT -q 30 ::: $OUT_DIR*/*.bam" >> $LOG.log
parallel -j $(echo `expr $N_THREADS \* 2`) --no-notice '(echo {}; infer_experiment.py -r ${GENOME_ANNOT_BED} -i {} -q 30)' ::: $(ls $OUT_DIR*/*[^e].bam) 2>&1 | tee -a $LOG.log

wait

# rsem_function () {
#     BAM_FILE=$1
#     BASE=$(basename -s ".toTranscriptome.bam" $BAM_FILE)
#     echo $BASE
#     mkdir -p ${OUT_DIR}${BASE}
#  
#     date > $OUT_DIR$BASE.tmp.log
#     pwd >> $OUT_DIR$BASE.tmp.log
#  
#     echo "----------------" >> $OUT_DIR$BASE.tmp.log
#     echo "Mapping $BASE to ref genome $GENOME_DIR_RSEM for RSEM" >> $OUT_DIR$BASE.tmp.log
#     (rsem-calculate-expression -p 2 \
#         --alignments \
#         --paired-end \
#         --strandedness reverse \
#         --estimate-rspd \
#         --output-genome-bam \
#         --sort-bam-by-coordinate \
#         --seed-length 23 \
#         $OUT_DIR$BASE/$BASE.toTranscriptome.bam \
#         $GENOME_DIR_RSEM \
#         $OUT_DIR$BASE/$BASE 2>&1 & echo $! >pid) | tee -a $OUT_DIR$BASE.tmp.log &
#     memory_script.sh $(cat pid) 5 &> $OUT_DIR$BASE.tmp.mem
#     wait
# }
# 
# export -f rsem_function
# parallel -j $N_THREADS --no-notice rsem_function ::: ${OUT_DIR}*/*.toTranscriptome.bam
# unset -f rsem_function
# 
# cat ${OUT_DIR}*.tmp.log > $LOG.log
# rm ${OUT_DIR}*.tmp.log
# 
# echo "-------------" > $LOG.mem
# echo "RSEM" >> $LOG.mem
# cat ${OUT_DIR}*.tmp.mem >> $LOG.mem
# rm ${OUT_DIR}*.tmp.mem
# wait

mkdir -p "${OUT_DIR}qualimap/rnaseq"
mkdir -p "${OUT_DIR}qualimap/bamqc"
# mkdir -p "${OUT_DIR}qualimap/multibamqc"
 
qualimap_rnaseq_function () {
    BAM_FILE=$1
    BASE=$(basename -s ".bam" $BAM_FILE)
    echo $BAM_FILE
    echo $BASE
    
    mkdir -p "${OUT_DIR}qualimap/rnaseq/$BASE"
    
    date > $OUT_DIR$BASE.tmp.log
    echo "-------------------" >> $OUT_DIR$BASE.tmp.log
    echo "Qualimap - rnaseq for $BASE" >> $OUT_DIR$BASE.tmp.log    
    # counts on uniquely mapped reads
    (JAVA_OPTS="-Djava.awt.headless=true -Djava.io.tmpdir=/home/diane/scratch" qualimap rnaseq -bam $BAM_FILE -gtf $GENOME_ANNOT_GTF -p strand-specific-reverse -pe -outdir "${OUT_DIR}qualimap/rnaseq/$BASE" --java-mem-size=4G 2>&1 & echo $! >pid) | tee -a $OUT_DIR$BASE.tmp.log
     
    # memory_script.sh $(cat pid) 5 &> $OUT_DIR$BASE.tmp.mem
    wait
}

export -f qualimap_rnaseq_function
parallel -j 4 --no-notice qualimap_rnaseq_function ::: $(ls ${OUT_DIR}*/*.bam)
unset -f qualimap_rnaseq_function

cat ${OUT_DIR}*.tmp.log >> $LOG.log
rm ${OUT_DIR}*.tmp.log
 
# echo "-------------" >> $LOG.mem
# echo "Qualimap Rnaseq" >> $LOG.mem
# cat ${OUT_DIR}*.tmp.mem >> $LOG.mem
# rm ${OUT_DIR}*.tmp.mem

wait

qualimap_bamqc_function () {    
    BAM_FILE=$1
    BASE=$(basename -s ".bam" $BAM_FILE)
    echo $BAM_FILE
    echo $BASE
    
    mkdir -p "${OUT_DIR}qualimap/bamqc/$BASE"
    
    date > $OUT_DIR$BASE.tmp.log
    echo "-------------------" >> $OUT_DIR$BASE.tmp.log
    echo "Qualimap - bamqc for $BASE" >> $OUT_DIR$BASE.tmp.log
    (JAVA_OPTS="-Djava.awt.headless=true -Djava.io.tmpdir=/home/diane/scratch" qualimap bamqc -bam $BAM_FILE -c -gd MOUSE -gff $GENOME_ANNOT_GTF -nt 4 -p strand-specific-reverse -outdir "${OUT_DIR}qualimap/bamqc/$BASE" --java-mem-size=4G 2>&1 & echo $! >pid)| tee -a $OUT_DIR$BASE.tmp.log &

    # memory_script.sh $(cat pid) 5 &> $OUT_DIR$BASE.tmp.mem 
    wait
    date >> $OUT_DIR$BASE.tmp.log
}

export -f qualimap_bamqc_function
parallel -j 6 --no-notice qualimap_bamqc_function ::: $(ls ${OUT_DIR}*/*.bam)
unset -f qualimap_bamqc_function

cat ${OUT_DIR}*.tmp.log >> $LOG.log
rm ${OUT_DIR}*.tmp.log

# echo "-------------" >> $LOG.mem
# echo "Qualimap Bamqc" >> $LOG.mem
# cat ${OUT_DIR}*.tmp.mem >> $LOG.mem
# rm ${OUT_DIR}*.tmp.mem

wait

# date >> $LOG.log
# echo "-------------------" >> $LOG.log
# echo "Qualimap - multibamqc" >> $LOG.log
# paste -d '\t' <(ls ${OUT_DIR}qualimap/bamqc) <(ls -d ${OUT_DIR}qualimap/bamqc/*) > data.tmp
# JAVA_OPTS="-Djava.awt.headless=true" qualimap multi-bamqc -d data.tmp -c -gff /opt/ngs_indexes/models/mm/mm10/gencode.vM6.refchrom.annotation.gtf -outdir "${OUT_DIR}qualimap/multibamqc/" --java-mem-size=4G 2>&1 | tee -a $LOG.log &

# wait
# rm data.tmp

date >> $LOG.log
echo "---------" >> $LOG.log
echo "Creating MultiQC reports" >> $LOG.log
mkdir -p ${OUT_DIR}multiqc
echo "multiqc $OUT_DIR  --outdir ${OUT_DIR}multiqc --interacive" >> $LOG.log
multiqc $OUT_DIR -f --outdir ${OUT_DIR}multiqc --interactive 2>&1 | tee -a $LOG.log

echo "Compressing all fastq files from $IN_DIR" >> $LOG.log
echo "parallel -j $(echo `expr $N_THREADS \* 2 `) --no-notice gzip -9 ::: ${IN_DIR}*${SUFF}" >> $LOG.log
parallel -j $(echo `expr $N_THREADS \* 2 `) --no-notice gzip -9 -f ::: ${IN_DIR}*${SUFF} 2>&1 | tee -a $LOG.log
 
date >> $LOG.log

wait 

unset IN_DIR OUT_DIR SUFF GENOME_DIR GENOME_ANNOT
rm pid
rm Log.out
rm Log.progress.out
rm -r _STARtmp
rm Aligned.out.sam

