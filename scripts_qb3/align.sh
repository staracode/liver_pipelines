#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o STAR_map_stdout                       #-- output directory (fill in)
#$ -e STAR_map_stderr                        #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
##$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=4G                  #-- submits on nodes with enough free memory (required)
#$ -l arch=linux-x64               #-- SGE resources (CPU type)
#$ -l netapp=10G,scratch=10G         #-- SGE resources (home and scratch disks)
#$ -l h_rt=24:00:00                #-- runtime limit (see above; this requests 24 hours)
#$ -t 1-15                        #-- remove first '#' to specify the number of

#READ in fastq file name and replace with trimmed fastq file name
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

#tasks=(0 /netapp/home/tfriedrich/Mattis/filtered//160715_I136_FCHCCTHBBXX_L6_WHHUMrkeRAADRAAPEI-209_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L1_HK500HUMyhuRAAORAAPEI-38_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L1_HK500HUMyhuRAAPRAAPEI-39_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L1_HK500HUMyhuRAAQRAAPEI-40_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L1_HK500HUMyhuRAARRAAPEI-41_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L1_HK500HUMyhuRAASRAAPEI-42_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L1_HK500HUMyhuRAATRAAPEI-43_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L1_HK500HUMyhuRAAWRAAPEI-1_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L2_HK500HUMyhuRABDRAAPEI-8_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L2_HK500HUMyhuRABERAAPEI-9_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L2_HK500HUMyhuRABFRAAPEI-10_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L2_HK500HUMyhuRABIRAAPEI-13_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L2_HK500HUMyhuRABJRAAPEI-14_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L2_HK500HUMyhuRABKRAAPEI-15_1_trimmed.fq.gz /netapp/home/tfriedrich/Mattis/filtered//170924_I89_CL100030081_L2_HK500HUMyhuRABLRAAPEI-16_1_trimmed.fq.gz)
fastq="${INPUT[$SGE_TASK_ID]}"
fastq_trim=`basename $fastq | sed 's/.fq.gz/_trimmed.fq.gz/'`
# Do I want to change fastq name to sample name here? 
prefix=`basename $fastq | sed 's/_trimmed.fq.gz/_algn/'`

# Output directory
if [ ! -d "$BIN_DIR" ]; then
	mkdir /netapp/home/tfriedrich/Mattis/algn/
fi 

echo $prefix
/netapp/home/tfriedrich/LiverCenter/software_source/STAR-2.6.0a/bin/Linux_x86_64/STAR --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFileNamePrefix   /netapp/home/tfriedrich/Mattis/algn/$prefix --readFilesCommand zcat --genomeDir /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/index/star  --readFilesIn $fastq_trim 