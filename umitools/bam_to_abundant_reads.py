#!/usr/bin/env python

import argparse
import sys
import pysam

global DEBUG
DEBUG = False


def read_to_key(read):
    read_chr_id = read.reference_id
    read5 = read.reference_start
    read_strand = ''
    if read.is_reverse:
        read_strand = '-'
    else:
        read_strand = '+'
    # We need abs() here as we do not know if r1 is the leftmost read or not
    tl = abs(read.template_length)
    k = ''.join((str(read_chr_id), ":", str(read5), "-", str(tl), read_strand))
    # print k
    return k


def get_reads_for_each_pos(fn):
    c = 0
    bam = pysam.AlignmentFile(fn, "rb")
    pos2read = {}
    for read in bam.fetch(reference="chrM"):
        c += 1
        if c % 100000 == 0:
            print >>sys.stderr, "Processed " + str(c) + " entries..."
        if read.is_read1:
            read_n = read.query_name
            # read_bc = read_n.split("_")[1]
            read_name = read_n.split("_")[0]
            k = read_to_key(read)
            if k in pos2read:
                pos2read[k].append(read_name)
            else:
                pos2read[k] = [read_name, ]
    return pos2read


def output_bam(fn, pos2read):
    all_read_names = {}
    for i in pos2read:
        for j in pos2read[i]:
            all_read_names[j] = 0
            
    bam = pysam.AlignmentFile(fn, "rb")
    out = pysam.AlignmentFile(fn + ".abundant_loci.bam", "wb", template=bam)
    c = 0
    for read in bam.fetch(reference="chrM"):
        c += 1
        if c % 100000 == 0:
            print >>sys.stderr, "Processed " + str(c) + " entries..."
        read_n = read.query_name
        # read_bc = read_n.split("_")[1]
        read_name = read_n.split("_")[0]
        if read_name in all_read_names:
            out.write(read)
    out.close()
            
    
def print_hist(pos2read):
    hist = {}
    for k in pos2read:
        l = len(pos2read[k])
        if l in hist:
            hist[l] += 1
        else:
            hist[l] = 1
    for i in range(max(hist.keys()) + 1):
        if i in hist:
            print str(i) + "\t" + str(hist[i])
        else:
            print str(i) + "\t" + str(0)

            
def main():
    parser = argparse.ArgumentParser(description='A pair of FASTQ files \
    are first reformatted using reformat_umi_fastq.py and then \
    is aligned to get the bam file. This script can parse the umi \
    barcode in the name of each read to mark duplicates.')
    parser.add_argument('-f', '--file', help='the input bam file',
                        required=False)
    parser.add_argument('-d', '--debug', help='turn on debug mode',
                        action="store_true")
    parser.add_argument('-o', '--output-bam', help='output the bam file containing reads of interest',
                        action="store_true")
    args = parser.parse_args()
    fn = args.file
    fn = "Zamore.RSQ.UMI.B6.brain.17dpp.fmt.x_rRNA.mm10g.sorted.bam"
    output_bam = False
    output_bam = args.output_bam
    bam_tmp = pysam.AlignmentFile(fn, "rb")


    if not bam_tmp.has_index():
        print >>sys.stderr, "Input file %s is being indexed." % fn
        pysam.index(fn)
    global DEBUG
    DEBUG = args.debug
    pos2read = get_reads_for_each_pos(fn)
    if output_bam:
        output_bam(fn, pos2read)
    print_hist(pos2read)


if __name__ == '__main__':
    main()
