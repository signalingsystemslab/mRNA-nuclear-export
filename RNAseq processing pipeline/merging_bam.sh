#!/bin/bash

EXP_FILE=$1
BASE=$(basename $EXP_FILE .txt)
export OUT_DIR="merged_bams_mm10_vM14_2/$BASE/"
LOG="merging_$BASE"

date 
mkdir -p ${OUT_DIR}
mkdir -p ${OUT_DIR}tmp
mkdir -p ${OUT_DIR}tmp_RG

for EXP in $(cat $EXP_FILE | cut -f 1 | sort | uniq)
do
    RUNS=$( cat $EXP_FILE | grep $EXP | cut -f2) 
#    FILES=$(cat $EXP_FILE | grep $EXP | cut -f2,3 --output-delimiter "_") 
    FILES=$(cat $EXP_FILE | grep $EXP | cut -f2,3 --output-delimiter "-") 

    # copy all raw bam files
    for i in $(seq 1 $(echo $FILES | wc -w))
    do
        RUN=$(echo $RUNS | cut -d" " -f $i)
        FILE=$(echo $FILES | cut -d" " -f $i)
#        cp /mnt/biggie/data2/diane/ribonomics/LPA_Naive_RNAseq/SxaQSEQs$RUN/06_filtering/b_deRibo/${FILE}.bam ${OUT_DIR}tmp
        cp lanes/SxaQSEQs$RUN/06_filtering/b_deRibo/${FILE}.bam ${OUT_DIR}tmp
        java -jar /opt/picard/picard.jar AddOrReplaceReadGroups I=${OUT_DIR}tmp/$FILE.bam O=${OUT_DIR}tmp_RG/$FILE.bam RGID=$FILE RGPL=ILLUMINA RGPU=$FILE RGSM=$EXP RGLB=$EXP
    
        rm ${OUT_DIR}tmp/$FILE.bam
    # CN center
    # DT date
    done

    # merge raw bam files to single experiment bam
    FILES=$(echo $FILES | sed "s@^@${OUT_DIR}tmp_RG/@g" | sed "s@ @ ${OUT_DIR}tmp_RG/@g" | sed 's@ @.bam @g' | sed 's@$@.bam@g')
   
    samtools merge -r ${OUT_DIR}$EXP.bam $FILES
    # samtools split -f '%!.%.' ${OUT_DIR}$EXP.bam    
    wait 
    
    rm ${OUT_DIR}tmp_RG/*.bam
done

# mkdir -p "${OUT_DIR}b_deRiboqualimap/rnaseq"
# mkdir -p "${OUT_DIR}b_deRibo/qualimap/bamqc"
# mkdir -p "${OUT_DIR}b_deRibo/qualimap/multibamqc"
# 
# qualimap_rnaseq_function () {
#     BAM_FILE=$1
#     BASE=$(basename -s ".bam" $BAM_FILE)
#     DIR=$(echo $BAM_FILE | cut -d"/" -f2)
#     # echo $BAM_FILE
#     # echo $BASE
#     # echo $DIR
# 
#     mkdir -p "${OUT_DIR}qualimap/rnaseq/$EXP"
# 
#     date >> $OUT_DIR$EXP.tmp.log
#     echo "-------------------" >> $OUT_DIR$EXP.tmp.log
#     echo "Qualimap - rnaseq for $EXP" >> $OUT_DIR$EXP.tmp.log
#     # counts on uniquely mapped reads
#     (JAVA_OPTS="-Djava.awt.headless=true" qualimap rnaseq -bam $BAM_FILE -gtf /opt/ngs_indexes/models/mm/mm10/gencode.vM6.refchrom.annotation.gtf -p strand-specific-reverse -pe -outdir "${OUT_DIR}qualimap/rnaseq/$EXP" --java-mem-size=4G 2>&1 & echo $! >pid) | tee -a $OUT_DIR$EXP.tmp.log &
# 
#     memory_script.sh $(cat pid) 5 &>> $OUT_DIR$EXP.tmp.mem
#     wait
# }
# 
# echo "rnaseq"
# export -f qualimap_rnaseq_function
# parallel -j 1 --no-notice qualimap_rnaseq_function ::: $(ls ${OUT_DIR}*[^{.junk},^{.ribo}].bam)
# unset -f qualimap_rnaseq_function
# 
# cat ${OUT_DIR}*/*.tmp.log >> $LOG.log
# rm ${OUT_DIR}*/*.tmp.log
# 
# echo "-------------" >> $LOG.mem
# echo "Qualimap Rnaseq" >> $LOG.mem
# cat ${OUT_DIR}*/*.tmp.mem >> $LOG.mem
# rm ${OUT_DIR}*/*.tmp.mem
# 
# wait
# 
# qualimap_bamqc_function () {
#     BAM_FILE=$1
#     EXP=$(basename -s ".bam" $BAM_FILE)
#     # echo $BAM_FILE
#     # echo $BASE
#     # echo $DIR
# 
#     mkdir -p "${OUT_DIR}/qualimap/bamqc/$EXP"
# 
#     date > $OUT_DIR$EXP.tmp.log
#     echo "-------------------" >> $OUT_DIR$EXP.tmp.log
#     echo "Qualimap - bamqc for $EXP" >> $OUT_DIR$EXP.tmp.log
#     (JAVA_OPTS="-Djava.awt.headless=true" qualimap bamqc -bam $BAM_FILE -c -gd MOUSE -gff /opt/ngs_indexes/models/mm/mm10/gencode.vM6.refchrom.annotation.gtf -nt 4 -p strand-specific-reverse -outdir "${OUT_DIR}qualimap/bamqc/$BASE" --java-mem-size=4G 2>&1 & echo $! >pid)| tee -a $OUT_DIR$EXP.tmp.log &
# 
#     memory_script.sh $(cat pid) 5 &> $OUT_DIR$EXP.tmp.mem
#     wait
#     date >> $OUT_DIR$EXP.tmp.log
# }
# 
# echo "bamqc"
# export -f qualimap_bamqc_function
# parallel -j 1 --no-notice qualimap_bamqc_function ::: $(ls ${OUT_DIR}${EXP}.bam)
# unset -f qualimap_bamqc_function
# 
# wait
# 
# date >> $LOG.log
# 
# echo "-------------------" >> $LOG.log
# echo "Qualimap - multibamqc" >> $LOG.log
# 
# paste -d '\t' <(ls ${OUT_DIR}/qualimap/bamqc) <(ls -d ${OUT_DIR}/qualimap/bamqc/*) > data.tmp
# JAVA_OPTS="-Djava.awt.headless=true" qualimap multi-bamqc -d data.tmp -c -gff /opt/ngs_indexes/models/mm/mm10/gencode.vM6.refchrom.annotation.gtf -outdir "${OUT_DIR}qualimap/multibamqc/" --java-mem-size=4G 2>&1 | tee -a $LOG.log &
# wait
# rm data.tmp
# 
# date >> $LOG.log
# echo "---------" >> $LOG.log
# echo "Creating MultiQC reports" >> $LOG.log
# 
# mkdir -p ${OUT_DIR}multiqc
# echo "multiqc $OUT_DIR  --outdir ${OUT_DIR}multiqc --interacive" >> $LOG.log
# multiqc ${OUT_DIR} -f --outdir ${OUT_DIR}multiqc --interactive 2>&1 | tee -a $LOG.log
# 
# date >> $LOG.log
# 
# wait

unset OUT_DIR
# rm pid

