#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o unmapped_fastq_stdout                       #-- output directory (fill in)
#$ -e unmapped_fastq_stderr                        #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
##$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=40G                  #-- submits on nodes with enough free memory (required)
#$ -l scratch=10G         #-- SGE resources (home and scratch disks)
#$ -l h_rt=24:00:00                #-- runtime limit (see above; this requests 24 hours)

DIR_NAME=$1

# File listing all the fastq file names and identifiers
FILE1=$2
samples=`wc -l $FILE1`
echo $samples
##$ -t 1-$samples  

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

# input file names 
fastq_file="${INPUT[$SGE_TASK_ID]}"
if [[ $fastq_file == *"fq.gz"* ]]; then
	fastq_unmapped_file=`echo $fastq_file | sed 's/.fq.gz/_unmapped.fq/'` 
	bam_file=`echo $fastq_file | sed 's/.fq.gz/_algn_Aligned.sortedByCoord.out.bam/'` 
else
	fastq_unmapped_file=`echo $fastq_file | sed 's/.fastq.gz/_unmapped.fq/'` 
	bam_file=`echo $fastq_file | sed 's/.fastq.gz/_algn_Aligned.sortedByCoord.out.bam/'` 
fi 

# Output directory
if [ ! -d "$BIN_DIR" ]; then
	mkdir /wynton/home/willenbring/tfriedrich/$DIR_NAME/unmapped/
fi

python LiverCenter/scripts_qb3/find_unmapped_reads.py $DIR_NAME/algn/$bam_file $DIR_NAME/fastq/$fastq_file  $DIR_NAME/unmapped/$fastq_unmapped_file
