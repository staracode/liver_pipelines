#$ -l mem_free=60G                  #-- submits on nodes with enough free memory (required)
#$ -l scratch=10G         #-- SGE resources (home and scratch disks)
#$ -l h_rt=24:00:00                #-- runtime limit (see above; this requests 24 hours)


localcores=16
localmem=50

#cd $TMPDIR
 
output_genome_name="mm10_with_markers_cellranger_index"
genome_fasta="LiverCenter/genomes/mm10/ucsc/mm10_with_markers.fa"
gtf_file="LiverCenter/genomes/mm10/ensembl/annotation/Mus_musculus.GRCm38.95_ucsc_formatted_rm_unusual_chr_with_markers.gtf"
cellranger-3.0.2/cellranger mkref \
	--genome="${output_genome_name}" \
	--fasta="${genome_fasta}" \
	--genes="${gtf_file}" \
	--nthreads="${localcores}" \
	--memgb="${localmem}"
