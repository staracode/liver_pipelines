files=`cut -d , -f 9 2019_05_02_Willenbring_rnaseq_tf_cholangiocyte_hepatocyte/GSE108315_SRA.csv`

echo $files

for fastq in $files; do 
	if [[ $fastq == *"SR"* ]]; then
		echo "LiverCenter/software_source/sratoolkit.2.9.6-centos_linux64/bin/fastq-dump " $fastq
	fi 
done
