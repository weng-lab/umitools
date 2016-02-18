#!/usr/bin/env python

# Usage: python umi_loci_with_duplicates.py -f test.sorted.bam 2>test.dup_marked.log

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
import sys
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
parser = argparse.ArgumentParser(description='A pair of FASTQ files are first reformatted using reformat_umi_fastq.py and then is aligned to get the bam file. This script can parse the umi barcode in the name of each read to mark duplicates.')
parser.add_argument('-f', '--file', help='the input bam file', required=True)
parser.add_argument('-p', '--processes', help='number of processes', required=False, type=int, default=8)
# parser.add_argument('-c', '--chromosome', help='the chromosome that you want to process', required=True)
# parser.add_argument('-a', '--add-tag', help='add a FM (five prime end of the mate) tag as the preprocessing step', action="store_true")
parser.add_argument('-d', '--debug', help='turn on debug mode', action="store_true")
parser.add_argument('-c', '--count', help='Count the number of raw reads for each locus (determined by pairs)', action="store_true")
# parser.add_argument('-l', '--read-length', help='if read length is given, it can be used to more accurately mark duplicates', type=int, default=-1)
# parser.add_argument('-o', '--output', help='the output file', required=True)
args = parser.parse_args()
processes = args.processes
global infile
infile = args.file
count_loc_flag = args.count
DEBUG = args.debug

# for kernprof
# kernprof -l -v ~/repo/tools/mark_duplicates_umi.py -f test300000.sorted.bam
# @profile
def find_loci_with_duplicates(infile, chromosome):
    bam = ps.AlignmentFile(infile, "rb")
    out = ps.AlignmentFile(infile + "." + chromosome + ".bam", "wb", template=bam)
    print >>sys.stderr, "Processing chromosome: " + chromosome
    r1fwd = -1
    r2fwd = -1
    # read5 + template length as the key
    counts_loc = {}
    out_counts_loc = open(infile + "." + chromosome + ".loc_count", "w")    
    # read5 + barcode + template length as the key
    counts_loc_umi = {}
    out_counts_loc_umi = open(infile + "." + chromosome + ".loc_umi_count", "w")
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
            read_id = str(read5) + "," + read_bc + "," + str(read.template_length) + ",r1"
            locus_id = str(read5) + "," + str(read.template_length) + ",r1"
            if locus_id in counts_loc:
                counts_loc[locus_id] += 1
            else:
                counts_loc[locus_id] = 1
            if read_id in counts_loc_umi:
                counts_loc_umi[read_id] += 1
            else:
                counts_loc_umi[read_id] = 1
            # if DEBUG:
            #     print >>sys.stderr, read_id + "\t" + str(read)
            # if read_id in r1fwd_ids:
            #     if DEBUG:
            #         print >>sys.stderr, "Found a duplicate (r1fwd): " + str(read)
            #     r1fwd_dup[read_n] = 'x'
            # else:
            #     r1fwd_ids[read_id] = 0
        elif not read.is_reverse and read.is_read2:
            # for a read that maps to the reverse strand, the read.template_length is still positive (not like those in the
            # bam file)
            read_id = str(read5) + "," + read_bc + "," + str(read.template_length) + ",r2"
            locus_id = str(read5) + "," + str(read.template_length) + ",r2"
            if locus_id in counts_loc:
                counts_loc[locus_id] +=1
            else:
                counts_loc[locus_id] = 1
            if read_id in counts_loc_umi:
                counts_loc_umi[read_id] += 1
            else:
                counts_loc_umi[read_id] = 1
            # if DEBUG:
            #     print >>sys.stderr, read_id + "\t" + str(read)
            # if read_id in r2fwd_ids:
            #     if DEBUG:
            #         print >>sys.stderr, read_n, "Found a duplicate (r2fwd)" + str(read)
            #     r2fwd_dup[read_n] = 'x'
            # else:
            #     r2fwd_ids[read_id] = 0
    bam.close()
    counts_loc_umi_cutoff = 1
    for i in counts_loc:
        print >> out_counts_loc, i + "\t" + str(counts_loc[i])
    for i in counts_loc_umi:
        if counts_loc_umi[i] >= counts_loc_umi_cutoff:
            print >>out_counts_loc_umi, i + "\t" + str(counts_loc_umi[i])
    out_counts_loc.close()
    out_counts_loc_umi.close()

def find_loci_with_duplicates_worker(chromosome):
    find_loci_with_duplicates(infile, chromosome)
    
if __name__ == "__main__":
    # mark_duplicates(bam, processes=processes)
    # mark_duplicates_helper(infile, processes=processes)
    bam_tmp = ps.AlignmentFile(infile, "rb")
    refs = bam_tmp.references
    p = Pool(processes)
    # Goddamned Pool does not accept lambda functions!
    p.map(find_loci_with_duplicates_worker, refs)
    # find_loci_with_duplicates(infile, "chr1")
    # mark_duplicates_worker("chr10")
    print("--- %s seconds ---" % (time.time() - start_time))
    
