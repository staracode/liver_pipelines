#wget ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz
#wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz
import os
import sys


# select genes on interest: 


cmd = "bin/bedtools2/bin/intersectBed -a Downloads/common_all_20180418\ \(1\).vcf  -b Downloads/Homo_sapiens.GRCh38.93.gff3"
import pprint
from BCBio.GFF import GFFExaminer
from BCBio import GFF

in_file = "/Users/tfriedrich/Downloads/Homo_sapiens.GRCh38.93.gff3"
#examiner = GFFExaminer()
#in_handle = open(in_file)
#pprint.pprint(examiner.parent_child_map(in_handle))
#in_handle.close()


#in_handle = open(in_file)
#for rec in GFF.parse(in_handle, target_lines=1000):
#    print rec
#in_handle.close()

genes = ["ACOT8", "NEF", "PEX5", "PTPN3", "IGHA1", "KLHL23", "PHOSPHO2-KLHL23", "PML", "PDHA1", "PTP4A1", "MTMR14", "NEDD1", "REL", "HCK", "ISG15", "MYC", "RNF2"]
genes_transcripts = [] # store ensembl transcript id here
filehandle = open(in_file, "r")
outhandle = open("/Users/tfriedrich/Downloads/Mattis_genes.gff3", "w")
for line in filehandle: 
    if line.startswith("#"): continue
    fields = line.rstrip("\n").split()
    #print fields
    annotation_type= fields[2]
    annotation = fields[8].split(";")
    if annotation_type == "mRNA": 
        geneName = annotation[2].split("-")[0].split("=")[1]
        if geneName in genes:
            #outhandle.write( line )
            transcript = annotation[0].replace("ID=transcript:", "")
            genes_transcripts.append(transcript)
    if annotation_type == "exon": 
        transcript = annotation[0].replace("Parent=transcript:", "")
        if transcript in genes_transcripts: 
            outhandle.write( line )
    
    cmd = "bin/bedtools2/bin/intersectBed -wao  -a   Downloads/Mattis_genes.gff3  -b Downloads/common_all_20180418\ \(1\).vcf"