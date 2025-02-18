#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o bowtie2_map_stdout                       #-- output directory (fill in)
#$ -e bowtie2_map_stderr                        #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
##$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=32G                  #-- submits on nodes with enough free memory (required)
#$ -l arch=linux-x64               #-- SGE resources (CPU type)
#$ -l netapp=10G,scratch=10G         #-- SGE resources (home and scratch disks)
#$ -l h_rt=24:00:00                #-- runtime limit (see above; this requests 24 hours)

# Directory where fastq file for project reside
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
	if [[ $fastq == *"fq.gz"* ]]; then
		fastq_trim2=`echo $fastq2 | sed 's/.fq.gz/_trimmed.fq.gz/'`
	else
		fastq_trim2=`echo $fastq2 | sed 's/.fastq.gz/_trimmed.fastq.gz/'`  
	fi 
fi 

# Do I want to change fastq name to sample name here? 
if [[ $fastq == *"fq.gz"* ]]; then
	prefix=`echo $fastq | sed 's/.fq.gz/_algn/'`
else
	prefix=`echo $fastq | sed 's/.fastq.gz/_algn/'`
fi

# Output directory
OUTPUTDIR=/netapp/home/tfriedrich/$DIR_NAME/algn/
if [ ! -d "$OUTPUTDIR" ]; then
	mkdir $OUTPUTDIR
fi 

FASTQ_TRIM_DIR=~/$DIR_NAME/filtered/

echo $OUTPUTDIR/$prefix
echo $prefix
echo $FASTQ_TRIM_DIR'/'$fastq_trim
echo $FASTQ_TRIM_DIR'/'$fastq_trim2 
INDEX=/netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/index/bowtie2/hg38_ucsc_2018_Nov_26

if [ "${#INPUT2[@]}" -gt 1 ]; then
	#bowtie2 --no-unal -x  /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/index/bowtie2/hg38_ucsc_2018_Nov_26   -1 $FASTQ_TRIM_DIR'/'$fastq_trim -2 $FASTQ_TRIM_DIR'/'$fastq_trim2  | samtools view -bS - > $OUTPUTDIR/$prefix
	bowtie2 --no-unal -x  /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/index/bowtie2/hg38_ucsc_2018_Nov_26   -1 $FASTQ_TRIM_DIR'/'$fastq_trim -2 $FASTQ_TRIM_DIR'/'$fastq_trim2  |  samtools view -bSu - | samtools sort - $OUTPUTDIR/$prefix

else
	#bowtie2 --no-unal -x  /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/index/bowtie2/hg38_ucsc_2018_Nov_26   -U $FASTQ_TRIM_DIR'/'$fastq_trim | samtools view -bS - > $OUTPUTDIR/$prefix
	bowtie2 --no-unal -x  /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/index/bowtie2/hg38_ucsc_2018_Nov_26   -U $FASTQ_TRIM_DIR'/'$fastq_trim | samtools view -bSu - | samtools sort - $OUTPUTDIR/$prefix
fi 
