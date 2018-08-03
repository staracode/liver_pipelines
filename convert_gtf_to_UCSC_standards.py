#converts ENSEMBL format to UCSC format for a GTF file
# Adds "chr" to each chromosome  
# Changes MT to M
import sys

all = map(str, range (0,50))
all.extend( ['MT', 'M', 'Y', 'X'] )
print (all)

gtf_file = "genomes/hg38/ensembl/annotation/Homo_sapiens.GRCh38.90.gtf"
track_chromosomes = {}
gtf_filehandle = open(gtf_file, 'r')
gtf_out = open(gtf_file.replace(".gtf", "_ucsc_formatted.gtf"), 'w')

for line in gtf_filehandle: 
    if line == "": gtf_out.write(line)
    fields = line.rstrip("\n").split("\t")
    if fields[0] in track_chromosomes: 
        inc = track_chromosomes[fields[0]]
        track_chromosomes[fields[0]] = inc + 1
    else: 
        track_chromosomes[fields[0]] = 1
    if fields[0] in all: 
        if fields[0] == "MT":
            chrom = "chrM"
        else:
            chrom = "chr" + fields[0]
        fields[0] = chrom
        new_line = "\t".join(fields)
        new_line = new_line + "\n"
        gtf_out.write(new_line) 
    else: 
        gtf_out.write(line)
for key in track_chromosomes:
    print"%s\t%s" % ( key,  track_chromosomes[key])