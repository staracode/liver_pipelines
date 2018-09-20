#wget ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz
#wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz
import os
import sys


# select genes on interest: 


cmd = "bin/bedtools2/bin/intersectBed -a Downloads/common_all_20180418\ \(1\).vcf  -b Downloads/Homo_sapiens.GRCh38.93.gff3"
import pprint
from BCBio.GFF import GFFExaminer

in_file = "/Users/tfriedrich/Downloads/Homo_sapiens.GRCh38.93.gff3"
examiner = GFFExaminer()
in_handle = open(in_file)
pprint.pprint(examiner.parent_child_map(in_handle))
in_handle.close()