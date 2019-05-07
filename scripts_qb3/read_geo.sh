#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o test_stdout                       #-- output directory (fill in)
#$ -e test_stderr                        #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
##$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=10G                  #-- submits on nodes with enough free memory (required)
#$ -l scratch=10G         #-- SGE resources (home and scratch disks)
#$ -l h_rt=24:00:00                #-- runtime limit (see above; this requests 24 hours)

# Directory where fastq file for project reside
DIR_NAME=$1

# File listing all the fastq file names and identifiers
FILE1=$2
#FILE1=2019_05_02_Willenbring_rnaseq_tf_cholangiocyte_hepatocyte/GSE108315_HW.csv

files=`cut -d , -f 9 ${FILE1}` 

for fastq in $files; do 
	if [[ $fastq == *"SR"* ]]; then
		echo "LiverCenter/software_source/sratoolkit.2.9.6-centos_linux64/bin/fastq-dump " $fastq " --outdir ${DIR_NAME}/fastq/ --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR_ID "  
	fi
done
