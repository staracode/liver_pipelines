#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o STAR_map_stdout                       #-- output directory (fill in)
#$ -e STAR_map_stderr                        #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
##$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=40G                  #-- submits on nodes with enough free memory (required)
#$ -l scratch=10G         #-- SGE resources (home and scratch disks)
#$ -l h_rt=24:00:00                #-- runtime limit (see above; this requests 24 hours)

# Directory where fastq file for project reside
DIR_NAME=$1

# File listing all the fastq file names and identifiers
FILE1=$2
samples=`wc -l $FILE1`
echo $samples
##$ -t 1-$samples  
SPECIES=$3  # hg38 or other (mm10)

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

# paired end 
fastq="${INPUT[$SGE_TASK_ID]}"
if [[ $fastq == *"fq.gz"* ]]; then
	fastq_trim=`echo $fastq | sed 's/.fq.gz/_trimmed.fq.gz/'` 
else
	fastq_trim=`echo $fastq | sed 's/.fastq.gz/_trimmed.fastq.gz/'` 
fi 

if [ "${#INPUT2[@]}" -gt 1 ]; then
	fastq2="${INPUT2[$SGE_TASK_ID]}"
	#fastq_trim2=`echo $fastq2 | sed 's/.fq.gz/_trimmed.fq.gz/'`  # TODO handle fastq.gz
	if [[ $fastq2 == *"fq.gz"* ]]; then
		fastq_trim2=`echo $fastq2 | sed 's/.fq.gz/_trimmed.fq.gz/'`  
	else
		fastq_trim2=`echo $fastq2 | sed 's/.fastq.gz/_trimmed.fastq.gz/'` 
	fi 
fi 
# Do I want to change fastq name to sample name here? 
if [[ $fastq == *"fq.gz"* ]]; then
	prefix=`echo $fastq | sed 's/.fq.gz/_algn_/'`
else
	prefix=`echo $fastq | sed 's/.fastq.gz/_algn_/'`
fi 

# Output directory
if [ ! -d "$BIN_DIR" ]; then
	mkdir /wynton/home/willenbring/tfriedrich/$DIR_NAME/algn/
fi 


FASTQ_TRIM_DIR=~/$DIR_NAME/filtered/

echo /wynton/home/willenbring/tfriedrich/$DIR_NAME/algn/$prefix
echo $prefix

if [ "${#INPUT2[@]}" -gt 1 ]; then
	if [[ $SPECIES == "hg38" ]]; then
		/wynton/home/willenbring/tfriedrich/LiverCenter/software_source/STAR-2.6.0a/bin/Linux_x86_64/STAR --outSAMtype BAM SortedByCoordinate  --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFileNamePrefix   /wynton/home/willenbring/tfriedrich/$DIR_NAME/algn/$prefix --readFilesCommand zcat --genomeDir /wynton/home/willenbring/tfriedrich/LiverCenter/genomes/hg38/ucsc/index/star  --readFilesIn $FASTQ_TRIM_DIR'/'$fastq_trim $FASTQ_TRIM_DIR'/'$fastq_trim2
	else
		/wynton/home/willenbring/tfriedrich/LiverCenter/software_source/STAR-2.6.0a/bin/Linux_x86_64/STAR --outSAMtype BAM SortedByCoordinate  --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFileNamePrefix   /wynton/home/willenbring/tfriedrich/$DIR_NAME/algn/$prefix --readFilesCommand zcat --genomeDir /wynton/home/willenbring/tfriedrich/LiverCenter/genomes/mm10/ucsc/index/star  --readFilesIn $FASTQ_TRIM_DIR'/'$fastq_trim $FASTQ_TRIM_DIR'/'$fastq_trim2
	fi 
else
	if [[ $SPECIES == "hg38" ]]; then
		/wynton/home/willenbring/tfriedrich/LiverCenter/software_source/STAR-2.6.0a/bin/Linux_x86_64/STAR --outSAMtype BAM SortedByCoordinate  --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFileNamePrefix   /wynton/home/willenbring/tfriedrich/$DIR_NAME/algn/$prefix --readFilesCommand zcat --genomeDir /wynton/home/willenbring/tfriedrich/LiverCenter/genomes/hg38/ucsc/index/star  --readFilesIn $FASTQ_TRIM_DIR'/'$fastq_trim 
	else
		/wynton/home/willenbring/tfriedrich/LiverCenter/software_source/STAR-2.6.0a/bin/Linux_x86_64/STAR --outSAMtype BAM SortedByCoordinate  --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFileNamePrefix   /wynton/home/willenbring/tfriedrich/$DIR_NAME/algn/$prefix --readFilesCommand zcat --genomeDir /wynton/home/willenbring/tfriedrich/LiverCenter/genomes/mm10/ucsc/index/star  --readFilesIn $FASTQ_TRIM_DIR'/'$fastq_trim 
	fi 
fi 
qstat -j $JOB_ID
