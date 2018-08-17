#move files
lftp --user tara.friedrich@ucsf.edu ftps://ftp.box.com

#check config

#quality metrics for fastq plus trimming
qsub -t 1-4 filter_fastq.sh config2.txt:

cd fastq
multiqc ./

cd filtered
multiqc ./



