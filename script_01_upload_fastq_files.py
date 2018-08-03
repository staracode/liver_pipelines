import pandas as pd
import os 

#read in fastq metadata table
#ls  /Volumes/My\ Passport/*/CleanData/*/*fq.gz > metadata_Aras_Mattis_2018_08_01.txt
metadata_file = "/Users/tarafriedrich/metadata_Aras_Mattis_2018_08_01.txt"
metadata = pd.read_csv(metadata_file, names=["location"])
#for this project only
full_filenames = metadata["location"].tolist()
filenames = [x.split("/")[-1] for x in full_filenames]
#print filenames
#os.system("mkdir ~/Mattis_fastq")
test = "rsync -avPh " + full_filenames[0].replace (" ", "\ ") + " tfriedrich@chef.compbio.ucsf.edu:~/Mattis_fastq"
print (test)
#how? pass2.compbio.ucsf.edu
#fastqc #quality control check 
#starr align to human 

#subread feature counts 

#downloaded hg38 on August 1, 2018 
#for x in {1..23} M X Y; do echo $x;   wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr'${x}'.fa.gz'; done;
#building index


#downloaded star rna-seq aligner version on August 1, 2018 
#wget https://github.com/alexdobin/STAR/archive/2.6.0a.tar.gz

# download ensembl hg38 genome and annotation on August 1, 2018
#rsync -avPh rsync://ftp.ensembl.org/ensembl/pub/current_embl/homo_sapiens  ./
#wget ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz



#copy bin
#source
star --runThreadN 8 --runMode genomeGenerate --genomeDir test/ --genomeFastaFiles /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr10.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr11.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr12.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr13.fa.gz,
/netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr14.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr15.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr16.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr17.fa.gz,
/netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr18.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr19.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr1.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr20.fa.gz,
/netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr21.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr22.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr2.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr3.fa.gz,
/netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr4.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr5.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr6.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr7.fa.gz,
/netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr8.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chr9.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chrM.fa.gz, /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chrX.fa.gz,
/netapp/home/tfriedrich/LiverCenter/genomes/hg38/ucsc/fasta/chrY.fa.gz --sjdbGTFfile /netapp/home/tfriedrich/LiverCenter/genomes/hg38/ensembl/annotation/Homo_sapiens.GRCh38.90_ucsc_formatted.gtf --sjdbOverhang 99