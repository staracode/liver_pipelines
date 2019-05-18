#taken from marton
from __future__ import division
import pysam
import multiprocessing
import datetime
import sys
import os
import argparse


# Command line argument parsing
parser = argparse.ArgumentParser(usage="python split_barcodes.py <options>", description='...')

parser.add_argument(
    "-b",
    default=None,
    dest='bam',
    action='store',
    help="Input bam file",
    required=True
)

parser.add_argument(
    "-c",
    default=None,
    dest='code',
    action='store',
    help="Barcode text file",
    required=True
)

parser.add_argument(
    "-n",
    default=None,
    dest='nproc',
    action='store',
    help="Number of processes",
    required=True
)

parser.add_argument(
    "-o",
    default=None,
    dest='output',
    action='store',
    help="Output directory",
    required=True
)
options = parser.parse_args(sys.argv[1:])

print ''
print datetime.datetime.now()

n_processes = int(options.nproc)


def split_list(v, n):
    ret = []
    for i in range(n):
        ret.append([])

    i = 0
    for x in v:
        ret[i].append(x)
        i += 1
        if i == n:
            i = 0

    return ret


bamfile = options.bam


chrs = ['chr{}'.format(x) for x in range(1, 20)] + ['chrX', 'chrY', 'chrM']
bam = pysam.Samfile(bamfile, "rb")
for r in bam.references:
    if r not in chrs:
        chrs.append(r)
chrom_sets = split_list(chrs, n_processes)

barcodes = []
for line in open(options.code):
    line = line.strip()
    if line == '':
        continue
    barcodes.append(line)


def process_bam(i):

    bam = pysam.Samfile(bamfile, "rb")
    for chrom in chrom_sets[i]:
        print 'Splitting: {}'.format(chrom)
        outfiles = {}
        for read in bam.fetch(chrom):
            read_barcode = read.get_tag('CR')
            if read_barcode in barcodes:
                if read_barcode not in outfiles:
                    outfiles[read_barcode] = pysam.AlignmentFile("{}/.{}_{}.bam".format(options.output, chrom, read_barcode), "wb", template=bam)
                outfiles[read_barcode].write(read)
        print 'Finished: {}'.format(chrom)


def merge_bams(i):

    for barcode in barcode_sets[i]:
        print 'Merging: {}'.format(barcode)
        os.system('samtools merge {}/.chr_{}.bam {}/.*_{}.bam'.format(options.output, barcode, options.output, barcode))
        os.system('samtools reheader {}/.header.sam {}/.chr_{}.bam > {}/{}.bam'.format(options.output, options.output, barcode, options.output, barcode))


processes = []
for i in range(n_processes):
    p = multiprocessing.Process(target=process_bam, args=(i, ))
    processes.append(p)

for i in range(n_processes):
    processes[i].start()

for i in range(n_processes):
    processes[i].join()

os.system('samtools view -H {} > {}/.header.sam'.format(bamfile, options.output))

barcode_sets = split_list(barcodes, n_processes)
processes = []
for i in range(n_processes):
    p = multiprocessing.Process(target=merge_bams, args=(i, ))
    processes.append(p)

for i in range(len(processes)):
    processes[i].start()

for i in range(len(processes)):
    processes[i].join()

os.system('rm {}/.chr*'.format(options.output))
os.system('rm {}/.header.sam'.format(options.output))

print datetime.datetime.now()
print ''
