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
#$ -t 1-15                        #-- remove first '#' to specify the number of

##tasks=(0 1bac 2xyz 3ijk 4abc 5def 6ghi 7jkl 8mno 9pqr 1stu )
##input="${tasks[$SGE_TASK_ID]}"

BIN_DIR=bin
FASTQ_DIR=~/Mattis/fastq/
FASTQ_TRIM_DIR=~/Mattis/filtered/
CMD=~/Mattis/cmds
if [ ! -d "$BIN_DIR" ]; then
	mkdir $BIN_DIR
fi
if [ ! -d "$FASTQ_TRIM_DIR" ]; then
	mkdir $FASTQ_TRIM_DIR
fi
if [ ! -d "$CMD" ]; then
	mkdir $CMD
fi 
cp -r ~/LiverCenter/software_source/FastQC $BIN_DIR
cp -r ~/LiverCenter/software_source/ExpressionAnalysis-ea-utils-bd148d4/clipper/fastq-mcf $BIN_DIR
cp ~/LiverCenter/pipeline_files/illumina.adapter.file.txt ./

#for fastq in  $FASTQ_DIR/*fq.gz; do 

tasks=(0 /netapp/home/tfriedrich/Mattis/fastq//160715_I136_FCHCCTHBBXX_L6_WHHUMrkeRAADRAAPEI-209_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L1_HK500HUMyhuRAAORAAPEI-38_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L1_HK500HUMyhuRAAPRAAPEI-39_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L1_HK500HUMyhuRAAQRAAPEI-40_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L1_HK500HUMyhuRAARRAAPEI-41_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L1_HK500HUMyhuRAASRAAPEI-42_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L1_HK500HUMyhuRAATRAAPEI-43_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L1_HK500HUMyhuRAAWRAAPEI-1_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L2_HK500HUMyhuRABDRAAPEI-8_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L2_HK500HUMyhuRABERAAPEI-9_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L2_HK500HUMyhuRABFRAAPEI-10_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L2_HK500HUMyhuRABIRAAPEI-13_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L2_HK500HUMyhuRABJRAAPEI-14_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L2_HK500HUMyhuRABKRAAPEI-15_1.fq.gz /netapp/home/tfriedrich/Mattis/fastq//170924_I89_CL100030081_L2_HK500HUMyhuRABLRAAPEI-16_1.fq.gz)
fastq="${tasks[$SGE_TASK_ID]}"
fastq_trim=`basename $fastq | sed 's/.fq.gz/_trimmed.fq.gz/'`
fastq_base=`basename $fastq | sed 's/.fq.gz//'`
	
ls 

echo $fastq_trim
#echo FastQC/fastqc $fastq > $CMD'/'$fastq_base'_cmd.sh'
#echo fastq-mcf illumina.adapter.file.txt $fastq -o $FASTQ_TRIM_DIR'/'$fastq_trim  >> $CMD'/'$fastq_base'_cmd.sh'
#echo FastQC/fastqc $fastq_trim >> $CMD'/'$fastq_base'_cmd.sh'

#/netapp/home/tfriedrich/LiverCenter/software_source/FastQC/fastqc $fastq 
#fastq-mcf illumina.adapter.file.txt $fastq -o $FASTQ_TRIM_DIR'/'$fastq_trim  
/netapp/home/tfriedrich/LiverCenter/software_source/FastQC/fastqc $FASTQ_TRIM_DIR'/'$fastq_trim 

#done

