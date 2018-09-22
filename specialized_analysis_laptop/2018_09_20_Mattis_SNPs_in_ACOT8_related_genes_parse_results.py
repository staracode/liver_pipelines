import sys

#parse results of intersection of VCF with regions of interest 
snps = "/Users/tfriedrich/Downloads/Mattis_genes_snps.txt"
outputfile = "/Users/tfriedrich/Downloads/Mattis_genes_snps_parsed.txt"
snps_filehandle = open(snps, "r")
outhandle = open(outputfile, "w")
for line in snps_filehandle: 
    if "ACOT8" in line:  # initially Aras Mattis was mainly interested in ACOT8 
        fields = line.rstrip("\n").split()
        chrom = "chr" + fields[0]
        snp_start = fields[10]
        snp_name = fields[11]
        ref = fields[12]
        alt = fields[13]
        anno = fields[16].split(";")  # I am having trouble parsing this
        outhandle.write("\t".join ([chrom, snp_start, snp_name,ref, alt, "\t".join(anno)]))
        outhandle.write("\n")
outhandle.close()
snps_filehandle.close()