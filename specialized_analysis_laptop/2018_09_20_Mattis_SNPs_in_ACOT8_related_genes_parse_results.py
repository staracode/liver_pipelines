import sys

#parse results of intersection of VCF with regions of interest 
snps = "/Users/tfriedrich/Documents/Mattis/Mattis_genes_snps.txt"
outputfile = "/Users/tfriedrich/Documents/Mattis/Mattis_genes_snps_parsed.txt"
outputfile2 = "/Users/tfriedrich/Documents/Mattis/Mattis_genes_snps_parsed_cytoscape.txt"
columns = ["chrom", "snp_position", "snp_id", "reference_allele", "alternative_allele", "gene_of_interest", "non-synonymous_stop_codon_present", "non-synonymous_missense_present", "non-sysnonymous_frameshit_present", "freqeuncies_major_minor_1000Genomes", "frequencies_major_minor_TOPMED"]
snps_filehandle = open(snps, "r")
outhandle = open(outputfile, "w")
outhandle2 = open(outputfile2, "w")
outhandle.write("\t".join(columns) + "\n")
for line in snps_filehandle: 
    if "ACOT8" in line:  # initially Aras Mattis was mainly interested in ACOT8 
        fields = line.rstrip("\n").split()
        chrom = "chr" + fields[0]
        snp_start = fields[10]
        snp_name = fields[11]
        ref = fields[12]
        alt = fields[13]
        [ freq1 , freq2 ] = [ "NA", "NA"]
        anno = fields[16].split(";")  # parse this line 
        for info in anno: 
            if "TOPMED" in info: 
                freq1 = info.replace("TOPMED=", "").split(",")
                
            if "CAF" in info: 
                freq2 = info.replace("CAF=", "").split(",")

            if "GENEINFO" in info: 
                gene = info.replace("GENEINFO=", "")
                genenames = gene.split('|')
                genes2 = [genename.split(":")[0] for genename in genenames]
                gene = "|".join(genes2)

        [NSN, NSM, NSF] = [0,0,0]
        if "NSN" in anno:
                NSN=1
            
        if "NSM" in anno:
                NSM=1
             
        if "NSF" in anno:
                NSF=1
         
            

        # column description 
        
        outhandle.write("\t".join ([chrom, snp_start, snp_name,ref, alt, "\t".join(map(str, [gene, NSN, NSM, NSF, freq1, freq2]))]))
        outhandle2.write("\t".join ([snp_name, "SNP", gene]))
        outhandle.write("\n")
        outhandle2.write("\n")
outhandle.close()
snps_filehandle.close()