#!/usr/bin/env python

import gzip
import sys
import argparse
from struct import unpack
import re

import umi

__author__ = "Yu Fu"
__license__ = "GPLv3"


def print2(a):
    print >>sys.stderr, a

    
def is_gzipped(filename):
    # 1F 8B 08 00 / gz magic number
    magic = ('\x1f', '\x8b', '\x08', '\x00')

    with open(filename, 'rb') as handle:
        s = unpack('cccc', handle.read(4))
        return s == magic

    
def process_sra_read(read, stats, ui):
    """This function processes one small RNA-seq read. Small RNA-seq
reads are always single end. Returns 4 empty strings if the read is dropped.
UMI information is passed to this function via the parameter ui
"""
    ret_name = ""
    ret_seq = ""
    ret_info = ""
    ret_qual = ""
    ret_bc = ""
    r_name = read.r_name
    r_seq = read.r_seq
    r_info = read.r_info
    r_qual = read.r_qual
    # mate = read.mate

    # print r_seq
    mm, bc = ui.check(r_seq)
    # print "---"

    if mm <= N_MISMATCH_ALLOWED_IN_UMI:
        ret_bc = bc
        ret_name = r_name
        ret_seq = ui.extract_insert(r_seq)
        ret_info = r_info
        ret_qual = ui.extract_insert_qual(r_qual)
    return ret_name, ret_seq, ret_info, ret_qual, ret_bc

    
def get_header_with_umi(header, umi):
    """This function inserts a UMI after the '@' symbol, making the downstream
    analysis easier
    """
    col = header.split(" ")
    return col[0] + "_" + umi + " " + col[1]


def is_good_phred(phred, qc):
    """Check if all the qualities are above qc, that is, if one
    nt has a quality score below qc, this function returns False
    """
    ret = True
    for i in phred:
        if i < qc:
            ret = False
            break
    return ret


# @profile
def main():
    desc = "A script to process reads in from UMI small RNA-seq. This script can handle \
    gzipped files transparently."
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input',
                        help='the input fastq file.', required=True)
    parser.add_argument('-o', '--output',
                        help='''the output fastq file containing reads that
                        are not duplicates''', required=True)
    parser.add_argument('-d', '--pcr-duplicate',
                        help='the output fastq file containing PCR duplicates',
                        required=True) 
    parser.add_argument('--reads-with-improper-umi', default="",
                        help='the output fastq file containing reads with improper UMIs. \
                        The default is to throw away these reads. This is for \
                        debugging purposes',
                        required=False)
   
    parser.add_argument('-v', '--verbose',
                        help='Also include detailed run info',
                        action="store_true")

    parser.add_argument('-5', '--umi-pattern-5',
                        help='''Set the UMI pattern at the 5\' end. Use ACGT for
                        fixed nt and N for variable nt in UMI. If there are
                        multiple patterns, separate them using comma''',
                        default='NNNCGANNNTACNNN,NNNATCNNNAGTNNN')
    
    parser.add_argument('-3', '--umi-pattern-3',
                        help='''Set the UMI pattern at the 3\' end. Use ACGT for
                        fixed nt and N for variable nt in UMI. If there are
                        multiple patterns, separate them using comma''',
                        default='NNNGTCNNNTAGNNN')
    # Quality cutoff is not necessary as the tools for trimming 3' small RNA-seq adapters
    # already do this
    # parser.add_argument('-q', '--quality',
    #                     help='Quality (phred quality score) cutoff for UMI. Default is 13, \
    #                     that is UMI with qualities >= 13 will be kept. This \
    #                     program assumes the phred quality scores in the fastq \
    #                     file are using sanger format',
    #                     required=False, type=int, default=13)
    global N_MISMATCH_ALLOWED_IN_UMI
    N_MISMATCH_ALLOWED_IN_UMI = 1
    c = 0
    args = parser.parse_args()
    
    fn = args.input
    verbose = args.verbose
    umi_pat5 = args.umi_pattern_5.split(",")
    umi_pat3 = args.umi_pattern_3.split(",")

    ui = umi.SraUmiInfo(umi_pat5, umi_pat3)

    print >>sys.stderr, "-" * 72
    print >>sys.stderr, "Summary of UMI patterns at 5' and 3' ends:"
    print >>sys.stderr, "-" * 72
    print >>sys.stderr, ui
    print >>sys.stderr, "-" * 72

    if re.search("\.gz|\.gzip", args.output):
        out = gzip.open(args.output, "wb", compresslevel=4)
    else:
        out = open(args.output, "w")
    if re.search("\.gz|\.gzip", args.pcr_duplicate):
        dup = gzip.open(args.pcr_duplicate, "wb", compresslevel=4)
    else:
        dup = open(args.pcr_duplicate, "w")
        
    if len(args.reads_with_improper_umi) != 0:
        if re.search("\.gz|\.gzip", args.reads_with_improper_umi):
            imp = gzip.open(args.reads_with_improper_umi, "wb",
                            compresslevel=4)
        else:
            imp = open(args.reads_with_improper_umi, "w")

    stats = umi.SraRunStats()
    if is_gzipped(fn):
        f = gzip.open(fn)
        print >>sys.stderr, "Input is gzipped."
    else:
        f = open(fn)

    ## Small RNA insert + umi should be unique; otherwise, it is a duplicate
    insert_umi = {}
    while True:
        c += 1
        if c % 4 == 1:
            r_name = f.readline().strip()
            if not r_name:
                break
        elif c % 4 == 2:
            r_seq = f.readline().strip()
        elif c % 4 == 3:
            r_info = f.readline().strip()
        else:
            r_qual = f.readline().strip()
            r = umi.SraRead(r_name, r_seq, r_info, r_qual, "r1")
            r_name_proc, r_seq_proc, r_info_proc, \
                r_qual_proc, r_bc = process_sra_read(r, stats, ui)
            if r_name_proc == "":
                stats["n_without_proper_umi"] += 1
                if len(args.reads_with_improper_umi) != 0:
                    print >>imp, r_name
                    print >>imp, r_seq
                    print >>imp, r_info
                    print >>imp, r_qual
            else:
                r_name_proc = get_header_with_umi(r_name_proc, r_bc)
                stats["n_with_proper_umi"] += 1
                k = r_seq + "_" + r_bc
                # print k
                if k in insert_umi:
                    stats["n_duplpicate"] += 1
                    print >>dup, r_name_proc
                    print >>dup, r_seq_proc
                    print >>dup, r_info_proc
                    print >>dup, r_qual_proc
                else:
                    stats["n_non_duplicate"] += 1
                    insert_umi[k] = True
                    print >>out, r_name_proc
                    print >>out, r_seq_proc
                    print >>out, r_info_proc
                    print >>out, r_qual_proc

    out.close()
    dup.close()
    
    if len(args.reads_with_improper_umi) != 0:
        imp.close()
    print >>sys.stderr, ""
    print >>sys.stderr, "Stats: "
    print >>sys.stderr, "Total input reads:\t" + str(c/4)
    print >>sys.stderr, "Reads dropped due to improper UMI:\t" + str(stats["n_without_proper_umi"])
    print >>sys.stderr, "Final proper read:\t" + str(stats["n_with_proper_umi"])
    print >>sys.stderr, "\tReads that are duplicates:\t" + str(stats["n_duplpicate"])
    print >>sys.stderr, "\tReads that are non-duplicates:\t" + str(stats["n_non_duplicate"])
        
    print2("")
    if verbose:
        pass

if __name__ == "__main__":
    main()
