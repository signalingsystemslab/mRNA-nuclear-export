#! /bin/bash

DB_DIR="/home/diane/software/meme/share/meme-5.0.3/db/motif_databases/RNA/"

mkdir -p 3utr_more_vs_less
mkdir -p 5utr_more_vs_less
mkdir -p 3utr_less_vs_more
mkdir -p 5utr_less_vs_more


scoring="totalhits"
evalue="1000"

for utr in "3utr" "5utr"
do
    hit_th=0.25
    method="spearman" 
    for direction in "more_vs_less" "less_vs_more"
       do 
       DIR="${utr}_${direction}/all_${method}_${scoring}_${hit_th}"
       if [ $direction == "more_vs_less" ]
       then
          ame --oc $DIR --method ${method}  --scoring ${scoring} --hit-lo-fraction $hit_th --evalue-report-threshold $evalue ${utr}/meme_${utr}_fasta_all.txt ${DB_DIR}Ray2013_rbp_Mus_musculus.meme ${DB_DIR}Ray2013_rbp_Homo_sapiens.meme
       else
          ame --oc $DIR --method ${method}  --scoring ${scoring} --hit-lo-fraction $hit_th --evalue-report-threshold $evalue ${utr}/meme_${utr}_fasta_all_reverse.txt ${DB_DIR}Ray2013_rbp_Mus_musculus.meme ${DB_DIR}Ray2013_rbp_Homo_sapiens.meme
       fi
    done
done
        

