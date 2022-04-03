## Naive

# Manual Check tracks and create text file of gene to export gene_list_pre_final_naive.txt
# Create gtf on the last 5kb
python3 get_5kb_gtf_with_dedup.py --gtf /opt/ngs_indexes/models/mm/mm10/gencode.vM14.annotation.gtf --out_gtf gencode.vM14_in_gene_list_final_naive_5kb.gtf --gene_list gene_list_final_naive.txt > Naive_5kb_remaining_bp.txt

# Add remaining bp to Excel file

# calculate count on all Naive bam files

N_THREADS=24

# take input from following directory
IN_DIR="merged_bams_mm10_vM14"
# write output files to following directory
OUT_DIR="merged_counts"
mkdir -p $OUT_DIR

# Genome annotation
GENOME_ANNOT="/opt/ngs_indexes/models/mm/mm10/gencode.vM14.annotation.gtf"

for DIR2 in $(ls -d ${IN_DIR}*/*Naive*/ | xargs -d '\n' -n 1 basename | sort | uniq )
do 
    echo $DIR2
    mkdir -p ${OUT_DIR}/${DIR2}"/D40000_LPAgenes_last5kb"
    
    featureCounts -T 8 -s 2 -S fr -p -B -C -P -d 40 -D 40000 -t exon -g gene_id -a /home/diane/annotations/gencode.vM14_in_gene_list_final_naive_5kb.gtf -o ${OUT_DIR}/${DIR2}/D40000_LPAgenes_last5kb/exon_counts.txt ${IN_DIR}*/${DIR2}/*.bam 2>&1 &
    
    wait
done

## LPA
# Manual Check tracks and create text file of gene to export gene_list_pre_final_naive.txt
# Create gtf on the last 5kb
python3 get_5kb_gtf_with_dedup.py --gtf /opt/ngs_indexes/models/mm/mm10/gencode.vM14.annotation.gtf --out_gtf gencode.vM14_in_gene_list_final_lpa_5kb.gtf --gene_list gene_list_final_lpa.txt > LPA_5kb_remaining_bp.txt

# calculate count on all LPA bam files

for DIR2 in $(ls -d ${IN_DIR}*/*LPA*/ | xargs -d '\n' -n 1 basename | sort | uniq )
do 
    echo $DIR2
    mkdir -p ${OUT_DIR}/${DIR2}"/D40000_LPAgenes_last5kb"
    
    featureCounts -T 8 -s 2 -S fr -p -B -C -P -d 40 -D 40000 -t exon -g gene_id -a /home/diane/annotations/gencode.vM14_in_gene_list_final_lpa_5kb.gtf -o ${OUT_DIR}/${DIR2}/D40000_LPAgenes_last5kb/exon_counts.txt ${IN_DIR}*/${DIR2}/*.bam 2>&1 &
    
    wait
done

