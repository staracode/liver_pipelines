#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o process_variant_rna_stdout                       #-- output directory (fill in)
#$ -e process_variant_rna_stderr                        #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
##$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=8G                  #-- submits on nodes with enough free memory (required)
#$ -l arch=linux-x64               #-- SGE resources (CPU type)
#$ -l scratch=10G         #-- SGE resources (home and scratch disks)
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

#variant calling rnaseq pipeline 
#https://software.broadinstitute.org/gatk/documentation/article.php?id=3891

# Do I want to change fastq name to sample name here? 
if [[ $fastq == *"fq.gz"* ]]; then
	suffix1=`echo $fastq | sed 's/.fq.gz/_algn.bam/'`
	suffix2=`echo $fastq | sed 's/.fq.gz/_algn_addRG.bam/'`
	suffix3=`echo $fastq | sed 's/.fq.gz/_algn_addRG_rmDup.bam/'`
	suffix4=`echo $fastq | sed 's/.fq.gz/_algn_addRG_rmDup.txt/'`
	suffix4_5=`echo $fastq | sed 's/.fq.gz/_algn_rmDup_addRG_SplitNCigar.bam/'`
	suffix5=`echo $fastq | sed 's/.fq.gz/_raw_snps_indels.g.vcf/'`
	suffix6=`echo $fastq | sed 's/.fq.gz/_raw_snps_indels_genotype.g.vcf/'`
	suffix7=`echo $fastq | sed 's/.fq.gz/_raw_snps_indels_genotype_filtered.g.vcf/'`
else
	suffix1=`echo $fastq | sed 's/.fastq.gz/_algn.bam/'`
	suffix2=`echo $fastq | sed 's/.fastq.gz/_algn_addRG.bam/'`
	suffix3=`echo $fastq | sed 's/.fastq.gz/_algn_addRG_rmDup.bam/'`
	suffix4=`echo $fastq | sed 's/.fastq.gz/_algn_addRG_rmDup.txt/'`
	suffix4_5=`echo $fastq | sed 's/.fastq.gz/_algn_rmDup_addRG_SplitNCigar.bam/'`
	suffix5=`echo $fastq | sed 's/.fastq.gz/_raw_snps_indels.g.vcf/'`
	suffix6=`echo $fastq | sed 's/.fastq.gz/_raw_snps_indels_genotype.g.vcf/'`
	suffix7=`echo $fastq | sed 's/.fastq.gz/_raw_snps_indels_genotype_filtered.g.vcf/'`
fi

# Output directory
OUTPUTDIR=/wynton/home/willenbring/tfriedrich/$DIR_NAME/algn/
if [ ! -d "$OUTPUTDIR" ]; then
	mkdir $OUTPUTDIR
fi 

echo $OUTPUTDIR/$prefix


#correctly naming samples is more important if I run the same sample multiple times based on this post:
#https://www.biostars.org/p/136302/

#descriptoin of reads gorups
#https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups

#extract ID and PU from the fastq file
# is my bam file sorted? 
# add read groups (so that I can merge later? )
bin/jdk1.8.0_191/bin/java -jar bin/picard.jar AddOrReplaceReadGroups \
       I=$OUTPUTDIR/$suffix1 \
       O=$OUTPUTDIR/$suffix2 \
       RGID=4 \
       RGLB=lib1 \
       RGPL=illumina \
       RGPU=unit1 \
       RGSM=20

#remove duplicates
bin/jdk1.8.0_191/bin/java -jar bin/picard.jar MarkDuplicates \
	#ASSUME_SORTED=true \
	I=$OUTPUTDIR/$suffix2 \ 
	O=$OUTPUTDIR/$suffix3 \
	M=$OUTPUTDIR/$suffix4 \
	CREATE_INDEX=true 

#samtools index $OUTPUTDIR/$suffix4

#java -jar picard.jar AddOrReplaceReadGroups I=star_output.sam O=rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample 
#java -jar picard.jar MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics 

bin/jdk1.8.0_191/bin/java -jar /wynton/home/willenbring/tfriedrich/bin/gatk-4.0.11.0/gatk-package-4.0.11.0-local.jar \
	GenomeAnalysisTK.jar \
	-T SplitNCigarReads \
	-R /wynton/home/willenbring/tfriedrich/LiverCenter/genomes/hg38/ucsc/all_chromosomes.fa \
	-I $OUTPUTDIR/$suffix3  \
	-o $OUTPUTDIR/$suffix4_5  \
	-rf ReassignOneMappingQuality \
	-RMQF 255 -RMQT 60 \
	-U ALLOW_N_CIGAR_READS

#haplotype caller
#https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.8.0/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php

bin/jdk1.8.0_191/bin/java -jar /wynton/home/willenbring/tfriedrich/bin/gatk-4.0.11.0/gatk-package-4.0.11.0-local.jar \
   HaplotypeCaller \
  -R /wynton/home/willenbring/tfriedrich/LiverCenter/genomes/hg38/ucsc/all_chromosomes.fa \
  -I $OUTPUTDIR/$suffix4_5 \
  -O $OUTPUTDIR/$suffix5 \
  --emit-ref-confidence GVCF

bin/jdk1.8.0_191/bin/java -jar /wynton/home/willenbring/tfriedrich/bin/gatk-4.0.11.0/gatk-package-4.0.11.0-local.jar \
   GenotypeGVCFs \
  -R /wynton/home/willenbring/tfriedrich/LiverCenter/genomes/hg38/ucsc/all_chromosomes.fa \
  -V $OUTPUTDIR/$suffix5 \
  -O $OUTPUTDIR/$suffix6

bin/jdk1.8.0_191/bin/java -jar /wynton/home/willenbring/tfriedrich/bin/gatk-4.0.11.0/gatk-package-4.0.11.0-local.jar \
 FilterVcf \
 --INPUT $OUTPUTDIR/$suffix6 \
 --OUTPUT $OUTPUTDIR/$suffix7 \
 --MIN_DP 30 --MIN_QD 30

