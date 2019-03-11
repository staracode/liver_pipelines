#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o filter_stdout                       #-- output directory (fill in)
#$ -e filter_stderr                        #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
##$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=4G                  #-- submits on nodes with enough free memory (required)
#$ -l scratch=10G         #-- SGE resources (home and scratch disks)
#$ -l h_rt=24:00:00                #-- runtime limit (see above; this requests 24 hours)

# Directory where fastq file for project reside
DIR_NAME=$1

# File listing all the fastq file names and identifiers
FILE1=$2
samples=`wc -l $FILE1`
echo $samples
##$ -t 1-$samples                       #-- remove first '#' to specify the number of

# Read in file containing fastq file locations. 
INPUT=(0)
INPUT2=(0)
while IFS=$'\t' read -r f1 f2 f3 #fastq, mate pair, label, replicate?
do
	printf 'FASTQ: %s, PAIR: %s\n' "$f1" "$f2" 
	# check whitespace
	# check that file is fastq.gz 
	# [ -z "$line" ] && continue
	INPUT+=($f1)
	if [[ $f2 == *"fq.gz"* ]]; then
		INPUT2+=($f2)  
	else
		if [[ $f2 == *"fastq.gz"* ]]; then 
			INPUT2+=($f2)
		else
			echo "unknown type of mate paired sequence file suffix"
		fi 
	fi
done <"$FILE1"
echo "${INPUT[@]}"

if [ "${#INPUT2[@]}" -gt 1 ]; then
	echo "mate"
	echo "${INPUT2[@]}"
fi

FASTQ_DIR=~/$DIR_NAME/fastq/
FASTQ_TRIM_DIR=~/$DIR_NAME/filtered/

if [ ! -d "$FASTQ_TRIM_DIR" ]; then
	mkdir $FASTQ_TRIM_DIR
fi

# Input Fastq Files
fastq="${INPUT[$SGE_TASK_ID]}"
if [[ $fastq == *"fq.gz"* ]]; then
	fastq_trim=`echo $fastq | sed 's/.fq.gz/_trimmed.fq.gz/'` 
else
	fastq_trim=`echo $fastq | sed 's/.fastq.gz/_trimmed.fastq.gz/'` 
fi 

if [ "${#INPUT2[@]}" -gt 1 ]; then
	fastq2="${INPUT2[$SGE_TASK_ID]}"
	if [[ $fastq2 == *"fq.gz"* ]]; then
		fastq_trim2=`echo $fastq2 | sed 's/.fq.gz/_trimmed.fq.gz/'`  
	else
		fastq_trim2=`echo $fastq2 | sed 's/.fastq.gz/_trimmed.fastq.gz/'` 
	fi 
fi

# Use this to debug issues
echo $PWD

echo $FASTQ_DIR'/'$fastq
echo $FASTQ_DIR'/'$fastq2

# Quality of data before trimming
if [ "${#INPUT2[@]}" -gt 1 ]; then
	/wynton/home/willenbring/tfriedrich/LiverCenter/software_source/FastQC/fastqc $FASTQ_DIR'/'$fastq $FASTQ_DIR'/'$fastq2
else
	/wynton/home/willenbring/tfriedrich/LiverCenter/software_source/FastQC/fastqc $FASTQ_DIR'/'$fastq
fi 

# Trim fastq filem
if [ "${#INPUT2[@]}" -gt 1 ]; then
	/wynton/home/willenbring/tfriedrich/LiverCenter/bin/bin/fastq-mcf ~/LiverCenter/pipeline_files/illumina.adapter.file.txt  $FASTQ_DIR'/'$fastq $FASTQ_DIR'/'$fastq2 -o $FASTQ_TRIM_DIR'/'$fastq_trim -o $FASTQ_TRIM_DIR'/'$fastq_trim2
else
	/wynton/home/willenbring/tfriedrich/LiverCenter/bin/bin/fastq-mcf ~/LiverCenter/pipeline_files/illumina.adapter.file.txt  $FASTQ_DIR'/'$fastq  -o $FASTQ_TRIM_DIR'/'$fastq_trim
fi 

# Quality of the data after trimming
if [ "${#INPUT2[@]}" -gt 1 ]; then
	/wynton/home/willenbring/tfriedrich/LiverCenter/software_source/FastQC/fastqc $FASTQ_TRIM_DIR'/'$fastq_trim $FASTQ_TRIM_DIR'/'$fastq_trim2
else
	/wynton/home/willenbring/tfriedrich/LiverCenter/software_source/FastQC/fastqc $FASTQ_TRIM_DIR'/'$fastq_trim
fi 

