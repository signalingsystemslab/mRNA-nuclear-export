#!/bin/bash
date

02_consolidated_michael.sh

if ls 02_consolidated/*.gz 2>/dev/null
then  
    parallel -j 24 --no-notice gzip -d ::: 02_consolidated/*.gz
fi
wait

02_fastqc.sh 02_consolidated/
wait

03_trimming_index_universal.sh 02_consolidated/  

wait

if ls 03_trimmed_index_universal/*.gz 2>/dev/null
then  
    parallel -j 24 --no-notice gzip -d -f ::: 03_trimmed_index_universal/*.gz
fi
wait

05_mapping.sh
wait

06_filtering.sh
wait

06_filtering_unique.sh
wait

07_tracks.sh
wait

07_tracks_unique.sh
wait

08_counts.sh
wait

# mkdir -p 09_multiqc
# multiqc . -d -f --outdir 09_multiqc --interactive

date
