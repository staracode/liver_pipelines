import gzip 
import sys
import pysam

bamfile = sys.argv[1]
#samtools view -h  chen_dataset/algn/NL_1_1_algn_Aligned.sortedByCoord.out.bam | head -1000 | samtools view -bS > test.bam
#bamfile = "test2.bam"
bam = pysam.AlignmentFile(bamfile, "rb")	
mapped_reads = [] #store aligned reads here
for line in bamfile:
	line =  str(line).rstrip("\n").split("\t")
	if line[0] in mapped_reads:
		continue
	else:
		mapped_reads.append(line[0])  #add read ID to list

fastqfile = sys.argv[2]
#zcat chen_dataset/fastq/NL_1_1.fq.gz | head -1000 | gzip > test.fq.gz
#fastqfile = "test.fq.gz"
fastq = gzip.open(fastqfile, "rb")

outfile = sys.argv[3]
out = open(outfile, "w")

nextLine=False
for line in fastq:
	if line.startswith("@"):
		nextLine = False
		if line.split(" ")[0].replace("@", "") in mapped_reads:
			continue
		else:
			#if this read is not in the bam file, print this line and the following three lines. 
			out.write( line )
			nextLine = True 
	else:
		if nextLine == True:
			out.write( line )
out.close()

