# mRNA-nuclear-export

## Step 1: RNAseq processing and QC
All files/scripts related to this step are listed under the `RNAseq processing pipeline` folder.  
RNAseq data were processed from qseq to counts using the pipeline `pipeline_ENCODE_ref.sh` ran on a Linux server ubuntu 16.04. The ouptut of this step are count files which are made available in the `Data > merged_counts` folder. We recommend starting after this step.  
For this step multiple scripts were used for each processing substeps:
* Conversion of qseq to fastq files see bash script`02_consolidated_michael.sh`
* Trimming adapter sequence using cudadapt, to see exact options used see `03_trimming_index_universal.sh`. Fastq files obtained after this step were deposit on ENCODE
* Reads were aligned using STAR software, to see exact options used see `05_trimming_index_universal.sh`
* Alignments were filtered to keep only reads mapped in proper pair and remove unmapped read/mates or failing vendor quality using samtools, to see exact options used see `06_filtering.sh`. 
* Alignment were filtered to only keep uniquely aligned reads using samtools. To see exact options used see `06_filtering_unique.sh`.
* Optional. Create bigwig track files to visualize sequencing with IGV using RSeQC `bamtowig.py` and USCS `wigToBigWig` tools, for exact options see `07_tracks.sh` and `07_tracks_unique.sh`
* Counts were generated using featureCounts. To see exact options used see `08_counts.sh`. Count files obtained from this step are available in the folder `Data > individual_counts`. 
* QC of the various steps was visualized using multiQC. Additionaly PCA using R was also done on all the technical replicates, see `QC > PCA.R`
* Bam files corresponding to technical replicates were merged into a single bam file using samtools, see `merging_bam.sh` for exact options. Bams files generated Bam files obtained after this step were deposited on ENCODE. Counts corresponding to merged bam files were generated using featureCounts. To see exact options used see `merged_counts.sh`. Count files obtained from this step are available in the folder `Data > merged_counts`. 

### Requirements for Step 1
Needed external softwares:
* python
* cutadapt
* STAR
* samtools
* RSeQC
* wigToBigWig
* multiQC
* R (ggplot2, edgeR)

## Step 2
DEG requires merged counts 
PCA DEGs and individual counts files
Create Datasets requires merged_counts i.e. after QC
Script to prepare dataset from raw count of step 1 is `Data > Create Datasets.R`. This will create R objects, with normalised counted using edgeR in the `Data` folder.
