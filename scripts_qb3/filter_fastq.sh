#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o filter_stdout                       #-- output directory (fill in)
#$ -e filter_stderr                        #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
##$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=4G                  #-- submits on nodes with enough free memory (required)
#$ -l arch=linux-x64               #-- SGE resources (CPU type)
#$ -l netapp=10G,scratch=10G         #-- SGE resources (home and scratch disks)
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
while IFS=: read -r f1 f2 
do
	printf 'Username: %s, Shell: %s\n' "$f1" "$f2" 
	# check whitespace
	# check that file is fastq.gz 
	# [ -z "$line" ] && continue
	INPUT+=($f1)
done <"$FILE1"
echo "${INPUT[@]}"

FASTQ_DIR=~/$DIR_NAME/fastq/
FASTQ_TRIM_DIR=~/$DIR_NAME/filtered/

if [ ! -d "$FASTQ_TRIM_DIR" ]; then
	mkdir $FASTQ_TRIM_DIR
fi

# Input Fastq Files
#tasks=(0 /netapp/home/tfriedrich/Mattis/fastq//160715_I136_FCHCCTHBBXX_L6_WHHUMrkeRAADRAAPEI-209_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L1_HK500HUMyhuRAAORAAPEI-38_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L1_HK500HUMyhuRAAPRAAPEI-39_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L1_HK500HUMyhuRAAQRAAPEI-40_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L1_HK500HUMyhuRAARRAAPEI-41_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L1_HK500HUMyhuRAASRAAPEI-42_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L1_HK500HUMyhuRAATRAAPEI-43_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L1_HK500HUMyhuRAAWRAAPEI-1_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L2_HK500HUMyhuRABDRAAPEI-8_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L2_HK500HUMyhuRABERAAPEI-9_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L2_HK500HUMyhuRABFRAAPEI-10_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L2_HK500HUMyhuRABIRAAPEI-13_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L2_HK500HUMyhuRABJRAAPEI-14_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L2_HK500HUMyhuRABKRAAPEI-15_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L2_HK500HUMyhuRABLRAAPEI-16_1.fq.gz)
fastq="${INPUT[$SGE_TASK_ID]}"
fastq_trim=`basename $fastq | sed 's/.fq.gz/_trimmed.fq.gz/'`
fastq_base=`basename $fastq | sed 's/.fq.gz//'`

# Use this to debug issues
echo $PWD
echo $fastq_trim

# Quality of data before trimming
/netapp/home/tfriedrich/LiverCenter/software_source/FastQC/fastqc $fastq 
# Trim fastq file
fastq-mcf ~/LiverCenter/pipeline_files/illumina.adapter.file.txt  $fastq -o $FASTQ_TRIM_DIR'/'$fastq_trim  
# Quality of the data after trimming
/netapp/home/tfriedrich/LiverCenter/software_source/FastQC/fastqc $FASTQ_TRIM_DIR'/'$fastq_trim 