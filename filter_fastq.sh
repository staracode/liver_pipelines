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
cp -r ~/LiverCenter/software_source/FastQC/fastqc $BIN_DIR
cp -r ~/LiverCenter/software_source/ExpressionAnalysis-ea-utils-bd148d4/clipper/fastq-mcf $BIN_DIR
cp ~/LiverCenter/pipeline_files/illumina.adapter.file.txt ./

#before trimming
for fastq in  $FASTQ_DIR/*fq.gz; do 

	fastq_trim=`basename $fastq | sed 's/.fq.gz/_trimmed.fq.gz/'`
	fastq_base=`basename $fastq | sed 's/.fq.gz//'`
	
	echo $fastq_trim
	echo fastqc $fastq > $CMD'/'$fastq_base'_cmd.sh'
	echo fastq-mcf illumina.adapter.file.txt $fastq -o $FASTQ_TRIM_DIR'/'$fastq_trim  >> $CMD'/'$fastq_base'_cmd.sh'
	echo fastqc $fastq_trim >> $CMD'/'$fastq_base'_cmd.sh'

done

