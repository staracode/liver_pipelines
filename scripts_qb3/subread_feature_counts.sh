#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o subread_fc_stdout                       #-- output directory (fill in)
#$ -e subread_fc_stderr                        #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
##$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=16G                  #-- submits on nodes with enough free memory (required)
#$ -l arch=linux-x64               #-- SGE resources (CPU type)
#$ -l netapp=10G,scratch=10G         #-- SGE resources (home and scratch disks)
#$ -l h_rt=24:00:00                #-- runtime limit (see above; this requests 24 hours)
##$ -t 1-15                        #-- remove first '#' to specify the number of


#tasks=(0 /netapp/home/tfriedrich/Mattis/filtered//160715_I136_FCHCCTHBBXX_L6_WHHUMrkeRAADRAAPEI-209_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L1_HK500HUMyhuRAAORAAPEI-38_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L1_HK500HUMyhuRAAPRAAPEI-39_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L1_HK500HUMyhuRAAQRAAPEI-40_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L1_HK500HUMyhuRAARRAAPEI-41_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L1_HK500HUMyhuRAASRAAPEI-42_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L1_HK500HUMyhuRAATRAAPEI-43_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L1_HK500HUMyhuRAAWRAAPEI-1_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L2_HK500HUMyhuRABDRAAPEI-8_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L2_HK500HUMyhuRABERAAPEI-9_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L2_HK500HUMyhuRABFRAAPEI-10_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L2_HK500HUMyhuRABIRAAPEI-13_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L2_HK500HUMyhuRABJRAAPEI-14_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L2_HK500HUMyhuRABKRAAPEI-15_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L2_HK500HUMyhuRABLRAAPEI-16_1_trimmed.fq.gz)
#fastq="${tasks[$SGE_TASK_ID]}"

if [ ! -d "$BIN_DIR" ]; then
	mkdir /netapp/home/tfriedrich/Mattis/counts
fi 

input_files=`ls Mattis/algn/*sam`
annotation_file=~/LiverCenter/genomes/hg38/ensembl/annotation/Homo_sapiens.GRCh38.90.gtf
output_file=~/Mattis/counts/readcounts.txt
~/LiverCenter/software_source/subread-1.6.2-source/bin/featureCounts -a $annotation_file -o $output_file $input_files
