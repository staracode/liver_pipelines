#!/bin/bash
 
localcores=16
localmem=60
 
# cellranger mkref --help
#output_genome_name="mm10_egfp_rabbit_v3"
output_genome_name="mm10_egfp_index"
genome_fasta="/proj_single-cell-3/genome/mm10_egfp.fa"
gtf_file="/proj_single-cell-3/genome/gtf_v3/ref-transcripts_egfp-filtered.gtf"
# is this filtered gtf file what victor meant when he said he didn't give me script that filtered something? 
/home/ubuntu/bin/cellranger mkref \
	--genome="${output_genome_name}" \
	--fasta="${genome_fasta}" \
	--genes="${gtf_file}" \
	--nthreads="${localcores}" \
	--memgb="${localmem}"




#!/bin/bash
illumina_bcl_dir="/proj_single-cell-3/data/runs/20190226-SATO/"
csv="early_cells_20190226-SATO.csv"
output_dir="/proj_single-cell-3/output_20190226-SATO"
localcores=16
localmem=64

bin/cellranger mkfastq \
    --run="${illumina_bcl_dir}" \
    --csv="${csv}" \
    --output-dir="${output_dir}" \
    --qc \ 
    --localcores="${localcores}" \
    --localmem="${localmem}"



#!/bin/bash 
 
localcores=16
localmem=60
 
fastqs="/proj_single-cell-3/output_20190325-SATO/HY3T7BCX2"
transcriptome="/proj_single-cell-3/mm10_egfp_index"
 
# Repeat this command per sample
for sample in KS_13  KS_15  KS_17
do
	/home/ubuntu/bin/cellranger count \
		--id="${sample}" \
		--sample="${sample}" \
		--fastqs="${fastqs}"  \
		--expect-cells=1000 \
		--transcriptome="${transcriptome}" \
		--localcores="${localcores}" \
		--localmem="${localmem}" \
		--nosecondary
done
