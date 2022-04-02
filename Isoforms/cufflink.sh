#!/bin/bash

# cufflinks on cytoplamic RNA

function cufflink_fun {
    NAME=$(basename --suffix .bam $1)
    DIRNAME=$( echo $NAME | sed 's@\.[0-9]\{3\}\.@\.@g')
    
    mkdir -p $DIRNAME
    mkdir -p $DIRNAME/$NAME
    cufflinks $1 -p 4 -G /opt/ngs_indexes/models/mm/mm10/gencode.vM14.annotation.gtf -o $DIRNAME/$NAME -u -b /opt/ngs_indexes/genomes/mm/mm10_no_alt_analysis_set_ENCODE.fasta  --library-type fr-firststrand 
}

export -f cufflink_fun
parallel -j 4 --no-notice cufflink_fun {} ::: ~/merged_bams_mm10_vM14/cytoplasmic.Naive.*/*.bam
parallel -j 4 --no-notice cufflink_fun {} ::: ~/merged_bams_mm10_vM14/cytoplasmic.LPA.*/*.bam
unset -f cufflink_fun
