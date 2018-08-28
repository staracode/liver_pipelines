#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o subread_fc_stdout                       #-- output directory (fill in)
#$ -e subread_fc_stderr                        #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
##$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=16G                  #-- submits on nodes with enough free memory (required)
#$ -l arch=linux-x64               #-- SGE resources (CPU type)
#$ -l netapp=10G,scratch=10G         #-- SGE resources (home and scratch disks)
#$ -l h_rt=24:00:00                #-- runtime limit (see above; this requests 24 hours)


# Directory where fastq file for project reside
DIR_NAME=$1

if [ ! -d "$BIN_DIR" ]; then
	mkdir /netapp/home/tfriedrich/$DIR_NAME/counts
fi 

input_files=`ls $DIR_NAME/algn/*bam`
annotation_file=~/LiverCenter/genomes/hg38/ensembl/annotation/Homo_sapiens.GRCh38.90.gtf
output_file=~/$DIR_NAME/counts/readcounts.txt
~/LiverCenter/software_source/subread-1.6.2-source/bin/featureCounts -p -O -M -T 1   -t exon -g gene_id -a $annotation_file -o $output_file $input_files
