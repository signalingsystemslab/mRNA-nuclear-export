#!/bin/bash

# make star genome index 
STAR --runThreadN 24 --runMode genomeGenerate \ 
--genomeDir /opt/ngs_indexes/star/mm10_no_alt_analysis_set_ENCODE.gencode_vM14_annotation.101bp \
--genomeFastaFiles /opt/ngs_indexes/genomes/mm/mm10_no_alt_analysis_set_ENCODE.fasta \
--sjdbGTFfile /opt/ngs_indexes/models/mm/mm10/gencode.vM14.annotation.gtf \ 
--sjdbOverhang 100 

