import pandas as pd
import os 

#read in fastq metadata table
#ls  /Volumes/My\ Passport/*/CleanData/*/*fq.gz > metadata_Aras_Mattis_2018_08_01.txt
metadata_file = "/Users/tarafriedrich/metadata_Aras_Mattis_2018_08_01.txt"
metadata = pd.read_csv(metadata_file, names=["location"])
#for this project only
full_filenames = metadata["location"].tolist()
filenames = [x.split("/")[-1] for x in full_filenames]
#os.system("mkdir ~/Mattis_fastq")
for name in full_filenames: 
    newname = name.replace(" ", "\ ")
    print "rsync -avPh " + newname + " tfriedrich@chef.compbio.ucsf.edu:~/Mattis_fastq"


#genomes/software_source/STAR-2.6.0a/bin/Linux_x86_64/STAR 
#subread feature counts 

#downloaded hg38 on August 1, 2018 ucsc
#for x in {1..23} M X Y; do echo $x;   wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr'${x}'.fa.gz'; done;
#building index

#downloaded star rna-seq aligner version on August 1, 2018 
#wget https://github.com/alexdobin/STAR/archive/2.6.0a.tar.gz

# download ensembl hg38 genome and annotation on August 1, 2018
#rsync -avPh rsync://ftp.ensembl.org/ensembl/pub/current_embl/homo_sapiens  ./
#wget ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz

#copy bin
#source
#genomes/software_source/STAR-2.6.0a/bin/Linux_x86_64/STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/index/star/ --genomeFastaFiles /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr10.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr11.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr12.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr13.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr14.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr15.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr16.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr17.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr18.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr19.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr1.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr20.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr21.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr22.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr2.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr3.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr4.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr5.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr6.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr7.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr8.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr9.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chrM.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chrX.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chrY.fa.gz --sjdbGTFfile /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ensembl/annotation/Homo_sapiens.GRCh38.90_ucsc_formatted.gtf --sjdbOverhang 48
#genomes/software_source/STAR-2.6.0a/bin/Linux_x86_64/STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/index/star/ --genomeFastaFiles /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr10.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr11.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr12.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr13.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr14.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr15.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr16.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr17.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr18.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr19.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr1.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr20.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr21.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr22.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr2.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr3.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr4.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr5.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr6.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr7.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr8.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr9.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chrM.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chrX.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chrY.fa.gz --sjdbGTFfile /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ensembl/annotation/Homo_sapiens.GRCh38.90_ucsc_formatted.gtf --sjdbOverhang 48


#fastqc #quality control check 
#./fastqc ~/Mattis_fastq/160715_I136_FCHCCTHBBXX_L6_WHHUMrkeRAADRAAPEI-209_1.fq.gz 
#starr align to human 
fastq-mcf 
multiqc 
alignment
subread feature counts
heatmap 
