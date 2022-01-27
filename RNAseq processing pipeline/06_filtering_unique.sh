#!/bin/bash
N_THREADS=8

# take input from following directory
export IN_DIR="06_filtering/b_deRibo/"
# write output files to foljlowing directory
export OUT_DIR="06_filtering/c_unique/"
#ribosomalRNA reference
export GENOME_ANNOT_GTF="/opt/ngs_indexes/models/mm/mm10/gencode.vM14.annotation.gtf"

mkdir -p $OUT_DIR

filter_function () { # for BASE in $(find $IN_DIR -maxdepth 1 -type d | cut -d"/" -f2)
# do
BAM_FILE=$1
BASE=$(basename -s ".bam" $BAM_FILE)

samtools view -b -h -F 2304 -@ 1 ${IN_DIR}${BASE}.bam -o ${OUT_DIR}${BASE}.bam &
wait

samtools index ${OUT_DIR}${BASE}.bam

wait

}

export -f filter_function
parallel -j $N_THREADS --no-notice filter_function ::: $(ls ${IN_DIR}*[^{.junk},^{.ribo}].bam)
unset -f filter_function

mkdir -p "${OUT_DIR}qualimap/rnaseq"
mkdir -p "${OUT_DIR}qualimap/bamqc"
# mkdir -p "${OUT_DIR}b_deRibo/qualimap/multibamqc"

qualimap_rnaseq_function () {
    BAM_FILE=$1
    BASE=$(basename -s ".bam" $BAM_FILE)
    # echo $BAM_FILE
    # echo $BASE
    # echo $DIR

    mkdir -p "${OUT_DIR}qualimap/rnaseq/$BASE"

    (JAVA_OPTS="-Djava.awt.headless=true -Djava.io.tmpdir=/home/diane/scratch" qualimap rnaseq -bam $BAM_FILE -gtf $GENOME_ANNOT_GTF -p strand-specific-reverse -pe -outdir "${OUT_DIR}qualimap/rnaseq/$BASE" --java-mem-size=4G 2>&1 ) 

    wait
}

export -f qualimap_rnaseq_function
parallel -j 4 --no-notice qualimap_rnaseq_function ::: $(ls ${OUT_DIR}*.bam)
unset -f qualimap_rnaseq_function

wait

qualimap_bamqc_function () {
    BAM_FILE=$1
    BASE=$(basename -s ".bam" $BAM_FILE)
    # echo $BAM_FILE
    # echo $BASE
    # echo $DIR

    mkdir -p "${OUT_DIR}qualimap/bamqc/$BASE"

    (JAVA_OPTS="-Djava.awt.headless=true -Djava.io.tmpdir=/home/diane/scratch" qualimap bamqc -bam $BAM_FILE -c -gd mm10 -gff $GENOME_ANNOT_GTF -nt 1 -p strand-specific-reverse -outdir "${OUT_DIR}qualimap/bamqc/$BASE" --java-mem-size=4G 2>&1 )

    wait
}

export -f qualimap_bamqc_function
parallel -j 4 --no-notice qualimap_bamqc_function ::: $(ls ${OUT_DIR}*.bam)
unset -f qualimap_bamqc_function

wait

mkdir -p ${OUT_DIR}multiqc
multiqc ${OUT_DIR} -f --outdir ${OUT_DIR}multiqc --interactive 2>&1 

wait

unset IN_DIR OUT_DIR 



