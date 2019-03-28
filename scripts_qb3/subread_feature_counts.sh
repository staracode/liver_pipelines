#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o subread_fc_stdout                       #-- output directory (fill in)
#$ -e subread_fc_stderr                        #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
##$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=32G                  #-- submits on nodes with enough free memory (required)
#$ -l scratch=10G         #-- SGE resources (home and scratch disks)
#$ -l h_rt=24:00:00                #-- runtime limit (see above; this requests 24 hours)


# Directory where fastq file for project reside
DIR_NAME=$1
SPECIES=$2

if [ ! -d "$BIN_DIR" ]; then
	mkdir /wynton/home/willenbring/tfriedrich/$DIR_NAME/counts
fi 


# TODO: throw an error if SPECIES isn't specified
input_files=`ls $DIR_NAME/algn/*bam`
if [[ $SPECIES == "hg38" ]]; then
	annotation_file=~/LiverCenter/genomes/hg38/ensembl/annotation/Homo_sapiens.GRCh38.90.gtf
else
	annotation_file=~/LiverCenter/genomes/mm10/ensembl/annotation/Mus_musculus.GRCm38.95_ucsc_formatted_rm_unusual_chr.gtf
fi 

output_file=~/$DIR_NAME/counts/readcounts.txt
~/LiverCenter/software_source/subread-1.6.2-source/bin/featureCounts -p -O -M -T 1   -t exon -g gene_id -a $annotation_file -o $output_file $input_files
