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
