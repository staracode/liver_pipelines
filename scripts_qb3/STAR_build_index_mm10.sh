#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o STAR_stdout                       #-- output directory (fill in)
#$ -e STAR_stderr                        #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
##$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=32G                  #-- submits on nodes with enough free memory (required)
#$ -l scratch=20G         #-- SGE resources (home and scratch disks)
#$ -l h_rt=24:00:00                #-- runtime limit (see above; this requests 24 hours)
##$ -t 8                        #-- remove first '#' to specify the number of
                                   #-- tasks if desired (see Tips section)
date
hostname

export TMPDIR=/scratch
export MYTMP=`mktemp -d`
cd $MYTMP

BINARY_DIR=/wynton/home/willenbring/tfriedrich/LiverCenter/software_source/STAR-2.6.0a/bin/Linux_x86_64/STAR
GTF_FILE=/wynton/home/willenbring/tfriedrich/LiverCenter/genomes/mm10/ensembl/annotation/Mus_musculus.GRCm38.95_ucsc_formatted_rm_unusual_chr.gtf
FASTA_DIR=/wynton/home/willenbring/tfriedrich/LiverCenter/genomes/mm10/ucsc/fasta
INDEX_DIR=star

cp -r $BINARY_DIR ./
cp -r $FASTA_DIR ./
cp $GTF_FILE ./

ls 

gunzip fasta/*
mkdir star 

./STAR --runMode genomeGenerate --genomeDir star/ --genomeFastaFiles fasta/chr10.fa fasta/chr11.fa fasta/chr12.fa fasta/chr13.fa fasta/chr14.fa fasta/chr15.fa fasta/chr16.fa fasta/chr17.fa fasta/chr18.fa fasta/chr19.fa fasta/chr1.fa fasta/chr2.fa fasta/chr3.fa fasta/chr4.fa fasta/chr5.fa fasta/chr6.fa fasta/chr7.fa fasta/chr8.fa fasta/chr9.fa fasta/chrM.fa fasta/chrX.fa fasta/chrY.fa --sjdbGTFfile Mus_musculus.GRCm38.95_ucsc_formatted_rm_unusual_chr.gtf --sjdbOverhang 100

cp -r $INDEX_DIR  /wynton/home/willenbring/tfriedrich/LiverCenter/genomes/mm10/ucsc/index/

rm -r fasta/

qstat -j $JOB_ID 
