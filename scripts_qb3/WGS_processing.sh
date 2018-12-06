#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o sort_bam_stdout                       #-- output directory (fill in)
#$ -e sort_bam_stderr                        #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
##$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=8G                  #-- submits on nodes with enough free memory (required)
#$ -l arch=linux-x64               #-- SGE resources (CPU type)
#$ -l netapp=10G,scratch=10G         #-- SGE resources (home and scratch disks)
#$ -l h_rt=14:00:00                #-- runtime limit (see above; this requests 24 hours)

# Directory where fastq file for project reside
DIR_NAME=$1

# File listing all the fastq file names and identifiers
FILE1=$2
samples=`wc -l $FILE1`
echo $samples
##$ -t 1-$samples  

# Read in file containing fastq file locations. 
INPUT=(0)
while IFS=$'\t' read -r f1 f2 f3 #fastq, mate pair, label, replicate?
do
	printf 'FASTQ: %s, PAIR: %s\n' "$f1" "$f2" 
	# check whitespace
	# check that file is fastq.gz 
	# [ -z "$line" ] && continue
	INPUT+=($f1)
done <"$FILE1"
echo "${INPUT[@]}"

# paired end 
fastq="${INPUT[$SGE_TASK_ID]}"

# Do I want to change fastq name to sample name here? 
if [[ $fastq == *"fq.gz"* ]]; then
	suffix1=`echo $fastq | sed 's/.fq.gz/_algn.bam/'`
	suffix2=`echo $fastq | sed 's/.fq.gz/_algn_rmDup.bam/'`
	suffix3=`echo $fastq | sed 's/.fq.gz/_algn_rmDup.txt/'`
	suffix4=`echo $fastq | sed 's/.fq.gz/_algn_rmDup_addRG.bam/'`
	suffix5=`echo $fastq | sed 's/.fq.gz/_raw_snps_indels.g.vcf/'`
else
	suffix1=`echo $fastq | sed 's/.fastq.gz/_algn.bam/'`
	suffix2=`echo $fastq | sed 's/.fastq.gz/_algn_rmDup.bam/'`
	suffix3=`echo $fastq | sed 's/.fastq.gz/_algn_rmDup.txt/'`
	suffix4=`echo $fastq | sed 's/.fastq.gz/_algn_rmDup_addRG.bam/'`
	suffix5=`echo $fastq | sed 's/.fastq.gz/_raw_snps_indels.g.vcf/'`
fi

# Output directory
OUTPUTDIR=/netapp/home/tfriedrich/$DIR_NAME/algn/
if [ ! -d "$OUTPUTDIR" ]; then
	mkdir $OUTPUTDIR
fi 

echo $OUTPUTDIR/$prefix

#remove duplicates
#bin/jdk1.8.0_191/bin/java -jar bin/picard.jar MarkDuplicates ASSUME_SORTED=true I=$OUTPUTDIR/$suffix1  O=$OUTPUTDIR/$suffix2  M=$OUTPUTDIR/$suffix3

# add read groups (so that I can merge later? )
#bin/jdk1.8.0_191/bin/java -jar bin/picard.jar AddOrReplaceReadGroups \
#       I=$OUTPUTDIR/$suffix2 \
#       O=$OUTPUTDIR/$suffix4 \
#       RGID=4 \
#       RGLB=lib1 \
#       RGPL=illumina \
#       RGPU=unit1 \
#       RGSM=20

#samtools index $OUTPUTDIR/$suffix4

bin/jdk1.8.0_191/bin/java -jar /netapp/home/tfriedrich/bin/gatk-4.0.11.0/gatk-package-4.0.11.0-local.jar \
   HaplotypeCaller \
  -R /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/all_chromosomes.fa \
  -I $OUTPUTDIR/$suffix4 \
  -O $OUTPUTDIR/$suffix5 \
  --emit-ref-confidence GVCF

