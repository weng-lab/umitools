#!/usr/bin/env python

# Usage: python mark_duplicates_umi.py -f test.sorted.bam >test.dup_marked.bam 2>test.dup_marked.log

__author__ = "Yu Fu"
__license__ = "GPLv3"

DEBUG = False
# A FASTQ file is first processed by reformat_umi_fastq.py to put the
# barcode info in the fasta header. Then it is aligned by aligner such as
# STAR. The bam file can then be processed by this script.
# It will go through the bam file and search for inserts that align to the same
# position and use UMI to mark all but one bam records as PCR duplicates
import pysam as ps
# from Bio.Seq import Seq
import sys
import argparse
import os
from multiprocessing import Process, Queue, Pool
import time

start_time = time.time()

# Two cases:
# ------------------------->
#                                     <-------------------------
#             r1                                   r2
#
# |
# r1fwd
#######################################################################
# ------------------------->
#                                     <-------------------------
#             r2                                   r1
# |
# r2fwd
######################################################################
# For a sorted bam file, in case 1, the 5' end of r1 and template length can determine the template;
# in case 2, the 5' end of r2 and template length can determine the template
# 

# for kernprof
# kernprof -l -v ~/repo/tools/mark_duplicates_umi.py -f test300000.sorted.bam
# @profile
def mark_duplicates(infile, chromosome):
    bam = ps.AlignmentFile(infile, "rb")
    out = ps.AlignmentFile(infile + "." + chromosome + ".bam", "wb", template=bam)
    print >>sys.stderr, "Processing chromosome: " + chromosome
    r1fwd = -1
    r2fwd = -1
    if count_loc_flag == True:
        # read5 + template length as the key
        counts_loc = {}
        # read5 + barcode + template length as the key
        out_counts_loc = open(infile + "." + chromosome + ".loc_count", "w")
    c = 0    
    # Unique read IDs: read5 + barcode + template length
    r1fwd_ids = {}
    r2fwd_ids = {}
    # Stores the reads to be marked as duplicates (only the 2nd and later ones are marked; the first one is not)
    r1fwd_dup = {}
    r2fwd_dup = {}
    for read in bam.fetch(reference=chromosome):
        c += 1
        if c % 100000 == 0:
            print >>sys.stderr, "Processed " + str(c) + " entries..."
        read_n = read.query_name
        read_bc = read_n.split("_")[1]
        read_chr = read.reference_id
        read5 = read.reference_start
        if DEBUG:
            print >>sys.stderr, "read_info: 5'end: " + str(read5)
            print >>sys.stderr, "read_info: Reference ID: " + str(read.reference_id)
            print >>sys.stderr, "read_info: Barcode: " + read_bc
            print >>sys.stderr, "read_info: Read: " + str(read)        
        if not read.is_reverse and read.is_read1:
            # We do not need to get abs(read.template_length) as it is guaranteed
            # to be positive here because we are only selecting the leftmost mate
            if count_loc_flag == True:
                locus_id = str(read5) + "," + str(read.template_length)
                if locus_id in counts_loc:
                    counts_loc[locus_id] += 1
                else:
                    counts_loc[locus_id] = 1
            read_id = (read5, read_bc, read.template_length)            
            if DEBUG:
                print >>sys.stderr, read_id + "\t" + str(read)
            if read_id in r1fwd_ids:
                if DEBUG:
                    print >>sys.stderr, "Found a duplicate (r1fwd): " + str(read)
                r1fwd_dup[read_n] = 'x'
            else:
                r1fwd_ids[read_id] = 0
        elif not read.is_reverse and read.is_read2:
            # for a read that maps to the reverse strand, the read.template_length is still positive (not like those in the
            # bam file)
            # read_id = str(read5) + read_bc + str(-read.template_length)
            # read_id = str(read5) + read_bc + str(read.template_length)
            # read_id = (read5, read_bc, read.template_length)
            # We do not need to get abs(read.template_length) as it is guaranteed
            # to be positive here because we are only selecting the leftmost mate
            read_id = str(read5) + read_bc + str(read.template_length)
            if count_loc_flag == True:
                locus_id = str(read5) + "," + str(read.template_length)
                if locus_id in counts_loc:
                    counts_loc[locus_id] +=1
                else:
                    counts_loc[locus_id] = 1
            if DEBUG:
                print >>sys.stderr, read_id + "\t" + str(read)
            if read_id in r2fwd_ids:
                if DEBUG:
                    print >>sys.stderr, read_n, "Found a duplicate (r2fwd)" + str(read)
                r2fwd_dup[read_n] = 'x'
            else:
                r2fwd_ids[read_id] = 0
    bam.close()
    # print >>sys.stderr, r1fwd_ids
    # print >>sys.stderr, r1fwd_dup
    # print >>sys.stderr, r2fwd_ids
    
    # The 2nd pass:
    bam = ps.AlignmentFile(infile, "rb")
    for read in bam.fetch(reference=chromosome):
        read_n = read.query_name
        if read_n in r1fwd_dup or read_n in r2fwd_dup:
            # Add 0x400 flags for PCR duplicates
            read.flag = read.flag | 0x400
        out.write(read)
    bam.close()
    out.close()
    if count_loc_flag:
        for locus_id in counts_loc:
            print >>out_counts_loc, chromosome + "\t" + locus_id + "\t" + str(counts_loc[locus_id])

def mark_duplicates_worker(chromosome):
    mark_duplicates(infile, chromosome)

def merge_bam(infile, refs):
    if infile.endswith(".bam"):
        prefix = infile[:-4]
    if prefix.endswith(".sorted"):
        prefix = prefix[:-7]
    # Open the input file to get the header template
    tmp = ps.AlignmentFile(infile, "rb")
    output_fn = prefix + ".deumi.sorted.bam"
    output = ps.AlignmentFile(output_fn, "wb", template=tmp)
    tmp.close()
    for chromosome in sorted(refs):
        fn = infile + "." + chromosome + ".bam"
        ps.index(fn)
        # print infile + "." + chromosome + ".bam"
        bam = ps.AlignmentFile(infile + "." + chromosome + ".bam", "rb")
        for read in bam.fetch(until_eof=True):
            output.write(read)
        bam.close()
        os.remove(fn)
        # Also remove the indices for the temp files
        os.remove(fn + ".bai")
    output.close()
    ps.index(output_fn)
    
def print2(a):
    print >>sys.stderr, a

def main():
    parser = argparse.ArgumentParser(description='A pair of FASTQ files are first reformatted using reformat_umi_fastq.py and then is aligned to get the bam file. This script can parse the umi barcode in the name of each read to mark duplicates.')
    parser.add_argument('-f', '--file', help='the input bam file', required=True)
    parser.add_argument('-p', '--processes', help='number of processes', required=False, type=int, default=8)
    # parser.add_argument('-c', '--chromosome', help='the chromosome that you want to process', required=True)
    # parser.add_argument('-a', '--add-tag', help='add a FM (five prime end of the mate) tag as the preprocessing step', action="store_true")
    parser.add_argument('-d', '--debug', help='turn on debug mode', action="store_true")
    parser.add_argument('-c', '--count', help='Count the number of raw reads for each locus (determined by pairs)', action="store_true", default=False)
    # parser.add_argument('-l', '--read-length', help='if read length is given, it can be used to more accurately mark duplicates', type=int, default=-1)
    # parser.add_argument('-o', '--output', help='the output file', required=True)
    args = parser.parse_args()
    processes = args.processes
    # True: output the sequences on the reference strand
    # False: output the sequences in the orignal direction
    global infile
    infile = args.file
    # bam = ps.AlignmentFile(infile, "rb")
    # readlength = args.read_length
    global count_loc_flag
    count_loc_flag = args.count
    DEBUG = args.debug

    # mark_duplicates(bam, processes=processes)
    # mark_duplicates_helper(infile, processes=processes)
    bam_tmp = ps.AlignmentFile(infile, "rb")
    if not bam_tmp.has_index():
        print >>sys.stderr, "Input file %s is being indexed." % infile
        ps.index(infile)
    refs = bam_tmp.references
    f = lambda x: mark_duplicates(infile, x)
    p = Pool(processes)
    # Goddamned Pool does not accept lambda functions!
    print2("Marking duplicates...")
    p.map(mark_duplicates_worker, refs)
    # mark_duplicates_worker("chr10")
    print2("Marking duplicates is done.")
    t1 = time.time()
    print2("--- %s seconds ---" % (t1 - start_time))
    print2("Merging from %d files..." % len(refs))
    merge_bam(infile, refs)
    print2("Merging is done")
    t2 = time.time()
    print2("--- %s seconds ---" % (t2-t1))

if __name__ == "__main__":
    main()
