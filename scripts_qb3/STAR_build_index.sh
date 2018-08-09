#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o STAR_stdout                       #-- output directory (fill in)
#$ -e STAR_stderr                        #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
##$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=32G                  #-- submits on nodes with enough free memory (required)
#$ -l arch=linux-x64               #-- SGE resources (CPU type)
#$ -l netapp=20G,scratch=20G         #-- SGE resources (home and scratch disks)
#$ -l h_rt=24:00:00                #-- runtime limit (see above; this requests 24 hours)
##$ -t 8                        #-- remove first '#' to specify the number of
                                   #-- tasks if desired (see Tips section)
date
hostname

export TMPDIR=/scratch
export MYTMP=`mktemp -d`
cd $MYTMP

BINARY_DIR=/netapp/home/tfriedrich/LiverCenter/bin
GTF_FILE=/netapp/home/tfriedrich/LiverCenter/genomes/hg38/ensembl/annotation/Homo_sapiens.GRCh38.90_ucsc_formatted.gtf
FASTA_DIR=/netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta
INDEX_DIR=star

cp -r $BINARY_DIR ./
cp -r $FASTA_DIR ./
cp $GTF_FILE ./

ls 

gunzip fasta/*
mkdir star 

#bin/STAR --runThreadN 8 --runMode genomeGenerate --genomeDir star/ --genomeFastaFiles fasta/chr10.fa fasta/chr11.fa fasta/chr12.fa fasta/chr13.fa fasta/chr14.fa fasta/chr15.fa fasta/chr16.fa fasta/chr17.fa fasta/chr18.fa fasta/chr19.fa fasta/chr1.fa fasta/chr20.fa fasta/chr21.fa fasta/chr22.fa fasta/chr2.fa fasta/chr3.fa fasta/chr4.fa fasta/chr5.fa fasta/chr6.fa fasta/chr7.fa fasta/chr8.fa fasta/chr9.fa fasta/chrM.fa fasta/chrX.fa fasta/chrY.fa --sjdbGTFfile Homo_sapiens.GRCh38.90_ucsc_formatted.gtf --sjdbOverhang 100
bin/STAR --runMode genomeGenerate --genomeDir star/ --genomeFastaFiles fasta/chr10.fa fasta/chr11.fa fasta/chr12.fa fasta/chr13.fa fasta/chr14.fa fasta/chr15.fa fasta/chr16.fa fasta/chr17.fa fasta/chr18.fa fasta/chr19.fa fasta/chr1.fa fasta/chr20.fa fasta/chr21.fa fasta/chr22.fa fasta/chr2.fa fasta/chr3.fa fasta/chr4.fa fasta/chr5.fa fasta/chr6.fa fasta/chr7.fa fasta/chr8.fa fasta/chr9.fa fasta/chrM.fa fasta/chrX.fa fasta/chrY.fa --sjdbGTFfile Homo_sapiens.GRCh38.90_ucsc_formatted.gtf --sjdbOverhang 100

cp -r $INDEX_DIR  /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/index/

rm -r bin
rm -r fasta/

qstat -j $JOB_ID 
