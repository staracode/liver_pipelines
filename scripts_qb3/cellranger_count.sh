#!/bin/bash 
#$ -l mem_free=60G                  #-- submits on nodes with enough free memory (required)
#$ -l scratch=10G         #-- SGE resources (home and scratch disks)
#$ -l h_rt=24:00:00                #-- runtime limit (see above; this requests 24 hours)
 
localcores=16
localmem=50
 
fastqs="FC_102418_fastqs"
transcriptome="LiverCenter/genomes/mm10_with_markers_cellranger_index/"
 
# Repeat this command per sample
for sample in FC_102418
do
	cellranger-3.0.2/cellranger count \
		--id="${sample}" \
		--sample="${sample}" \
		--fastqs="${fastqs}"  \
		--expect-cells=7000 \
		--transcriptome="${transcriptome}" \
		--localcores="${localcores}" \
		--localmem="${localmem}" \
		--nosecondary 
		#--chemistry="SC3Pv2"
done
