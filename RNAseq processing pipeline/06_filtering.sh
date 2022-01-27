#!/bin/bash

#  06_filtering.sh
#  
#
#  Created by Diane on 10/26/16.
#

# discard non-mappers, chrM and chrY reads, improperly paired reads and reads mapping to rRNA coordinates

N_THREADS=24

# take input from following directory
export IN_DIR="05_alignment/"
# write output files to foljlowing directory
export OUT_DIR="06_filtering/"
#ribosomalRNA reference
export RRNA_ANNOT="/home/diane/mm10_rRNA.bed"
export GENOME_ANNOT_GTF="/opt/ngs_indexes/models/mm/mm10/gencode.vM14.annotation.gtf"

LOG="06_filtering"

# create OUT_DIR if it doesn't exist
mkdir -p $OUT_DIR
mkdir -p ${OUT_DIR}a_filter/
mkdir -p ${OUT_DIR}b_deRibo/
 
pwd > $LOG.log
date >> $LOG.log

filter_function () { # for BASE in $(find $IN_DIR -maxdepth 1 -type d | cut -d"/" -f2)
# do
BAM_FILE=$1
BASE=$(basename -s ".bam" $BAM_FILE)
echo "" >> $OUT_DIR$BASE.tmp.log 
echo "-----------------" >> $OUT_DIR$BASE.tmp.log
date >> $OUT_DIR$BASE.tmp.log
echo "Filtering $BASE:" >> $OUT_DIR$BASE.tmp.log
echo "Reads or mates unmapped: $(samtools view -f 12 -@ 1 ${IN_DIR}${BASE}/${BASE}.bam | cut -f1 | sort | uniq | wc -l)" >> $OUT_DIR$BASE.tmp.log # 12  
echo "Reads fail platform/vendor quality checks: $(samtools view -f 512 -@ 1 ${IN_DIR}${BASE}/${BASE}.bam | cut -f1 | sort | uniq | wc -l)" >> $OUT_DIR$BASE.tmp.log # 512  
echo "Reads on chrX: $(samtools view -@ 1 ${IN_DIR}${BASE}/${BASE}.bam chrX chrX_GL456233_random | cut -f1 | sort | uniq | wc -l)" >> $OUT_DIR$BASE.tmp.log
echo "Reads on chrY: $(samtools view -@ 1 ${IN_DIR}${BASE}/${BASE}.bam chrY chrY_JH584301_random chrY_JH584300_random chrY_JH584303_random chrY_JH584302_random | cut -f1 | sort | uniq | wc -l)" >> $OUT_DIR$BASE.tmp.log
echo "Reads on chrM: $(samtools view -@ 1 ${IN_DIR}${BASE}/${BASE}.bam chrM | cut -f1 | sort | uniq | wc -l)" >> $OUT_DIR$BASE.tmp.log
#echo "Properly mapped pairs: $(samtools view -f 2 -@ 1 ${IN_DIR}${BASE}/${BASE}.bam | cut -f1 | sort | uniq | wc -l)" >> $OUT_DIR$BASE.tmp.log # 2  

wait 

# Removing read unmapped, mate unmapped, read fails platform/vendor quality checks -F 524; keep only properly mapped pairs -f 2;  
# awk: remove chrM and chrY
# # # echo "cat $(samtools view -H ${IN_DIR}${BASE}/${BASE}.bam) $(samtools view -F 2828 -f 2 -q 30 -@ 2 ${IN_DIR}${BASE}/${BASE}.bam | awk '{if($3 != "chrM" && $3 != "chrY"){print $0}}') | samtools view -@ 1 -b - > ${OUT_DIR}a_filter/${BASE}.bam"
samtools view -H ${IN_DIR}${BASE}/${BASE}.bam > ${OUT_DIR}a_filter/${BASE}_header.sam
samtools view -F 524 -f 2 -@ 1 ${IN_DIR}${BASE}/${BASE}.bam | awk '{if($3 != "chrM" && !match($3,"chrY")){print $0}}' | cat ${OUT_DIR}a_filter/${BASE}_header.sam - |  samtools view -@ 1 -b - > ${OUT_DIR}a_filter/${BASE}.bam & 

# memory_script.sh $! 1 > $OUT_DIR$BASE.tmp.mem 
wait 

rm ${OUT_DIR}a_filter/${BASE}_header.sam

# echo "Pairs Remaining: $(samtools view -@ 1 ${OUT_DIR}a_filter/${BASE}.bam | cut -f1 | sort | uniq | wc -l)" >> $OUT_DIR$BASE.tmp.log # 2  
# echo "Pairs Remaining (unique): $(samtools view -q 255 -@ 1 ${OUT_DIR}a_filter/${BASE}.bam | cut -f1 | sort | uniq | wc -l)" >> $OUT_DIR$BASE.tmp.log # 2 

# wait 

echo "Indexing BAM files" >> $OUT_DIR$BASE.tmp.log
samtools index ${OUT_DIR}a_filter/${BASE}.bam

wait 

echo "Filter ribosomal RNA for $BASE:" >> $OUT_DIR$BASE.tmp.log
echo "split_bam.py -i ${OUT_DIR}a_filter/${BASE}.bam -r $RRNA_ANNOT -o ${OUT_DIR}b_deRibo/${BASE}" >> $OUT_DIR$BASE.tmp.log

(split_bam.py -i ${OUT_DIR}a_filter/${BASE}.bam -r $RRNA_ANNOT -o ${OUT_DIR}b_deRibo/${BASE} 2>&1 & echo $! >pid) | tee -a $OUT_DIR$BASE.tmp.log &
# memory_script.sh $(cat pid) 1 >> $OUT_DIR$BASE.tmp.mem

wait 
    
mv ${OUT_DIR}b_deRibo/${BASE}.ex.bam ${OUT_DIR}b_deRibo/${BASE}.bam
mv ${OUT_DIR}b_deRibo/${BASE}.in.bam ${OUT_DIR}b_deRibo/${BASE}.ribo.bam

echo "Pairs in rRNA: $(samtools view -@ 1 ${OUT_DIR}b_deRibo/${BASE}.ribo.bam | cut -f1 | sort | uniq | wc -l)" >> $OUT_DIR$BASE.tmp.log # 2  
echo "Pairs in remaining after filtering: $(samtools view -@ 1 ${OUT_DIR}b_deRibo/${BASE}.bam | cut -f1 | sort | uniq | wc -l)" >> $OUT_DIR$BASE.tmp.log # 2
echo "Pairs in remaining after filtering (unique): $(samtools view -q 255 -@ 1 ${OUT_DIR}b_deRibo/${BASE}.bam | cut -f1 | sort | uniq | wc -l)" >> $OUT_DIR$BASE.tmp.log # 2 
  
wait

echo "Indexing BAM files" >> $OUT_DIR$BASE.tmp.log
samtools index ${OUT_DIR}b_deRibo/${BASE}.bam

wait 

# Potential tRNA
echo "Potential tRNA: $(samtools view -@ 1 -L /home/diane/mm10_tRNA.bed -c ${OUT_DIR}b_deRibo/${BASE}.bam)"

date >> $OUT_DIR$BASE.tmp.log

#done
}

export -f filter_function
parallel -j $N_THREADS --no-notice filter_function ::: $(ls ${IN_DIR}*/*.bam)
unset -f filter_function

cat ${OUT_DIR}*.tmp.log >> $LOG.log
rm ${OUT_DIR}*.tmp.log

# echo "-------------" >> $LOG.mem
# echo "Filtering" >> $LOG.mem
# cat ${OUT_DIR}*.tmp.mem >> $LOG.mem
# rm ${OUT_DIR}*.tmp.mem

# mkdir -p "${OUT_DIR}a_filter/qualimap/rnaseq"
# mkdir -p "${OUT_DIR}a_filter/qualimap/bamqc"
# mkdir -p "${OUT_DIR}a_filter/qualimap/multibamqc"

mkdir -p "${OUT_DIR}b_deRibo/qualimap/rnaseq"
mkdir -p "${OUT_DIR}b_deRibo/qualimap/bamqc"
# mkdir -p "${OUT_DIR}b_deRibo/qualimap/multibamqc"

qualimap_rnaseq_function () {
    BAM_FILE=$1
    BASE=$(basename -s ".bam" $BAM_FILE)
    DIR=$(echo $BAM_FILE | cut -d"/" -f2)
    # echo $BAM_FILE
    # echo $BASE
    # echo $DIR
    
    mkdir -p "${OUT_DIR}${DIR}/qualimap/rnaseq/$BASE"

    date >> $OUT_DIR$DIR/$BASE.tmp.log
    echo "-------------------" >> $OUT_DIR$DIR/$BASE.tmp.log
    echo "Qualimap - rnaseq for ${DIR}/$BASE" >> $OUT_DIR$DIR/$BASE.tmp.log
    # counts on uniquely mapped reads
    (JAVA_OPTS="-Djava.awt.headless=true -Djava.io.tmpdir=/home/diane/scratch" qualimap rnaseq -bam $BAM_FILE -gtf $GENOME_ANNOT_GTF -p strand-specific-reverse -pe -outdir "${OUT_DIR}${DIR}/qualimap/rnaseq/$BASE" --java-mem-size=4G 2>&1 & echo $! >pid) | tee -a $OUT_DIR$DIR/$BASE.tmp.log &

#    memory_script.sh $(cat pid) 5 &>> $OUT_DIR$DIR/$BASE.tmp.mem
    wait
}
echo "rnaseq"
export -f qualimap_rnaseq_function
parallel -j 6 --no-notice qualimap_rnaseq_function ::: $(ls ${OUT_DIR}b_deRibo/*[^{.junk},^{.ribo}].bam)
unset -f qualimap_rnaseq_function

cat ${OUT_DIR}*/*.tmp.log >> $LOG.log
rm ${OUT_DIR}*/*.tmp.log

# echo "-------------" >> $LOG.mem
# echo "Qualimap Rnaseq" >> $LOG.mem
# cat ${OUT_DIR}*/*.tmp.mem >> $LOG.mem
# rm ${OUT_DIR}*/*.tmp.mem

wait

qualimap_bamqc_function () {
    BAM_FILE=$1
    BASE=$(basename -s ".bam" $BAM_FILE)
    DIR=$(echo $BAM_FILE | cut -d"/" -f2)
    # echo $BAM_FILE
    # echo $BASE
    # echo $DIR

    mkdir -p "${OUT_DIR}${DIR}/qualimap/bamqc/$BASE"

    date > $OUT_DIR$DIR/$BASE.tmp.log
    echo "-------------------" >> $OUT_DIR$DIR/$BASE.tmp.log
    echo "Qualimap - bamqc for $DIR/$BASE" >> $OUT_DIR$DIR/$BASE.tmp.log
    (JAVA_OPTS="-Djava.awt.headless=true -Djava.io.tmpdir=/home/diane/scratch" qualimap bamqc -bam $BAM_FILE -c -gd MOUSE -gff $GENOME_ANNOT_GTF -nt 4 -p strand-specific-reverse -outdir "${OUT_DIR}${DIR}/qualimap/bamqc/$BASE" --java-mem-size=4G 2>&1 & echo $! >pid)| tee -a $OUT_DIR$DIR/$BASE.tmp.log &

#   memory_script.sh $(cat pid) 5 &> $OUT_DIR$DIR/$BASE.tmp.mem
    wait
    date >> $OUT_DIR$DIR/$BASE.tmp.log
}

echo "bamqc"
export -f qualimap_bamqc_function
parallel -j 6 --no-notice qualimap_bamqc_function ::: $(ls ${OUT_DIR}b_deRibo/*[^{.junk},^{.ribo}].bam)
unset -f qualimap_bamqc_function

# echo "-------------" >> $LOG.mem
# echo "Qualimap Bamqc" >> $LOG.mem
# cat ${OUT_DIR}*/*.tmp.mem >> $LOG.mem
# rm ${OUT_DIR}*/*.tmp.mem

cat ${OUT_DIR}*/*.tmp.log >> $LOG.log
rm ${OUT_DIR}*/*.tmp.log

wait

# date >> $LOG.log
# echo "-------------------" >> $LOG.log
# echo "Qualimap - multibamqc" >> $LOG.log
# # paste -d '\t' <(ls ${OUT_DIR}a_filter/qualimap/bamqc) <(ls -d ${OUT_DIR}a_filter/qualimap/bamqc/*) > data.tmp
# # JAVA_OPTS="-Djava.awt.headless=true" qualimap multi-bamqc -d data.tmp -c -gff /opt/ngs_indexes/models/mm/mm10/gencode.vM6.refchrom.annotation.gtf -outdir "${OUT_DIR}a_filter/qualimap/multibamqc/" --java-mem-size=4G 2>&1 | tee -a $LOG.log &

# # wait

# paste -d '\t' <(ls ${OUT_DIR}b_deRibo/qualimap/bamqc) <(ls -d ${OUT_DIR}b_deRibo/qualimap/bamqc/*) > data.tmp
# JAVA_OPTS="-Djava.awt.headless=true" qualimap multi-bamqc -d data.tmp -c -gff /opt/ngs_indexes/models/mm/mm10/gencode.vM6.refchrom.annotation.gtf -outdir "${OUT_DIR}b_deRibo/qualimap/multibamqc/" --java-mem-size=4G 2>&1 | tee -a $LOG.log &
# wait
# rm data.tmp

date >> $LOG.log
echo "---------" >> $LOG.log
echo "Creating MultiQC reports" >> $LOG.log
# mkdir -p ${OUT_DIR}a_filter/multiqc
# echo "multiqc ${OUT_DIR}a_filter/ --outdir ${OUT_DIR}a_filter/multiqc --interacive" >> $LOG.log
# multiqc ${OUT_DIR}a_filter/ -f --outdir ${OUT_DIR}a_filter/multiqc --interactive 2>&1 | tee -a $LOG.log

mkdir -p ${OUT_DIR}b_deRibo/multiqc
echo "multiqc $OUT_DIR  --outdir ${OUT_DIR}b_deRibo/multiqc --interacive" >> $LOG.log
multiqc ${OUT_DIR}b_deRibo/ -f --outdir ${OUT_DIR}b_deRibo/multiqc --interactive 2>&1 | tee -a $LOG.log

date >> $LOG.log

wait

unset IN_DIR OUT_DIR RRNA_ANNOT
rm pid
