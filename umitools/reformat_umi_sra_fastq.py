#!/usr/bin/env python

__author__ = "Yu Fu"
__license__ = "GPLv3"

import operator
import gzip
import sys
import argparse
from struct import unpack


def print2(a):
    print >>sys.stderr, a


class SraRunStats():
    """It records the info for run (for small RNA-seq data)
"""
    def __init__(self):
        self.stats = {}
        self.stats["all_umi_locator"] = {}
        self.stats["n_without_locator"] = 0
        self.stats["n_with_proper_umi"] = 0
        self.stats["n_without_proper_umi"] = 0
        # Those reads with N's before GGG
        self.stats["n_with_ambiguous_umi"] = 0
        # Those with A/C/G after GGG
        self.stats["n_with_wrong_padding"] = 0
        # Those having bad quality UMI w/ GGG and T and w/o N's
        self.stats["n_bad_quality_umi"] = 0
        # The rest
        self.stats["n_good_reads"] = 0
        # This stores counts of all UMIs for good reads only
        self.stats["good_umi_locator"] = {}
        # padding (usually T, or it can be A,T,C,G
        self.stats["padding"] = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
        # For ligation bias
        self.stats["ligation"] = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
    def __getitem__(self, key):
        return self.stats[key]
    def __setitem__(self, key, value):
        self.stats[key] = value
    def summary(self):
        print2("-" * 80)
        print2("Total:\t" + str((c-1)/4) )
        print2("Reads w/o locator:\t" + str(self["n_without_locator"]) )
        print2("Reads w/ locator:\t" + str(self["n_with_locator"]) )
        print2("-" * 80 )
        print2("Reads w/ N's in UMI:\t" + str(self["n_with_ambiguous_umi"]) )
        print2("Reads w/ wrong padding nt:\t" + str(self["n_with_wrong_padding"]) )
        print2("Reads w/ low-quality UMI:\t" + str(self["n_bad_quality_umi"]) )
        print2("-" * 80)
        print2("Reads w/ proper UMI:\t" + str(self["n_good_reads"]) )
        print2("-" * 80)
    def umi_padding_usage(self):
        print2("UMI padding usage for good reads")
        for i in self.stats["good_umi_locator"]:
            print2( i + " " + str(self.stats["good_umi_locator"][i]) )
        print2("UMI padding usage for all reads")
        my_sorted = sorted(self.stats["all_umi_locator"].items(), key=operator.itemgetter(1), reverse=True)
        for i in my_sorted:
            print2( str(i[0]) + " " + str(i[1]) )
    def padding_usage(self):
        for i in sorted(self.stats["padding"]):
            print2(i + "\t" + str(self.stats["padding"][i]))
    def ligation_bias(self):
        for i in sorted(self.stats["ligation"]):
            print2(i + "\t" + str(self.stats["ligation"][i]))
        
        
class Read():
    """ Informtion for one read including r_name, r_seq, r_info, r_qual, mate and so on...
"""
    def __init__(self, r_name, r_seq, r_info, r_qual, mate):
        self.r_name = r_name
        self.r_seq = r_seq
        self.r_info = r_info
        self.r_qual = r_qual
        self.mate = mate
        # self.my_umi_locator = r_seq[umi_len : umi_len + umi_locator_len]
        # self.my_umi = my_umi = r_seq[0:umi_len]
        # self.my_umi_qual = r_qual[0:umi_len]


def is_gzipped(filename):
    # 1F 8B 08 00 / gz magic number
    magic = ('\x1f', '\x8b', '\x08', '\x00')

    with open(filename, 'rb') as handle:
        s = unpack('cccc', handle.read(4))
        return s==magic


# @profile
def process_sra_read(read, stats):
    """This function processes one small RNA-seq read. Small RNA-seq
reads are always single end. Returns 4 empty strings if the read is dropped
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
    mate = read.mate
    if DEBUG:
        print '-' * 80
        print mate
        print r_name
        print "Original qual:\t" + r_qual
        print "Original read:\t" + r_seq
        print "UMI:\t\t" + r_seq[0: umi_len]
    has_proper_umi_5 = True
    has_proper_umi_3 = True
    for i in range(0, len(umi_pat_5)):
        # When the SRA read does not have proper umi
        if umi_pat_5[i] != "N":
            if r_seq[i] != umi_pat_5[i]:
                has_proper_umi_5 = False
                print "Read dropped due to improper 5\' UMI"
                print "Correct nt %s; nt in read %s" % (umi_pat_5[i], r_seq[i])
                print r_seq
                break
        else:
            ret_bc += r_seq[i]
    # for i in range(0, len(umi_pat_3)):
    #     ii = len(r_seq) - len(umi_pat_3)
    #     if umi_pat_3[i] != "N":
    #         if r_seq[ii] != umi_pat_3[i]:
    #             has_proper_umi_3 = False
    #             print "Read dropped due to improper 3\' UMI"
    #             print r_seq
    #             break
    #     else:
    #         ret_bc += r_seq[ii]
            
    if has_proper_umi_5 and has_proper_umi_3:
        ret_name = "TODO " + r_name
        ret_seq = "TODO" + r_seq
        ret_info = "TODO " + r_info
        ret_qual = "TODO" + r_qual
        # my_umi_qual = read.my_umi_qual # r_qual[0:umi_len]
        # if my_umi.find("N") == -1:
        #     padding = r_seq[umi_len + umi_locator_len]
        #     if padding in umi_downstream:
        #         stats[mate]["padding"][padding] += 1
        #         if is_good_phred(my_umi_qual, qc):
        #             my_seq = r_seq[umi_len + umi_locator_len + umi_downstream_len: ]
        #             stats[mate]["n_good_reads"] += 1
        #             stats[mate]["ligation"][my_seq[0]] += 1
        #             # I do not modify read name here,
        #             # since we need to store barcodes from both reads and put the
        #             # concatenated barcode into both reads to make downstream analysis easier
        #             # ret_name = get_header_with_umi(r_name, my_umi)
        #             ret_name = r_name
        #             ret_seq =  my_seq
        #             ret_info = r_info
        #             ret_qual = r_qual[umi_len + umi_locator_len + umi_downstream_len: ]
        #             ret_bc = my_umi
        #             if my_umi_locator in stats[mate]["good_umi_locator"]:
        #                 stats[mate]["good_umi_locator"][my_umi_locator] += 1
        #             else:
        #                 stats[mate]["good_umi_locator"][my_umi_locator] = 1
        #         else:
        #             stats[mate]["n_bad_quality_umi"] += 1
        #             if DEBUG:
        #                 print "*Found one with low quality UMI:\t" + r_seq
        #     else:
        #         stats[mate]["n_with_wrong_padding"] += 1
        #         if DEBUG:
        #             print "*Found one with wrong padding nucleotide (A/C/G):\t" + r_seq
        # else:
        #     if DEBUG == True:
        #         print "*Found one with N's in UMI:\t" + r_seq
        #     stats[mate]["n_with_ambiguous_umi"] += 1
        return ret_name, ret_seq, ret_info, ret_qual, ret_bc
    else:
        if DEBUG:
            print "*Found one without proper umi:\t" + r_seq
        stats["n_without_proper_umi"] += 1
        # return ret_name, ret_seq, ret_info, ret_qual, ret_bc
        return ["", ] * 5

    
def get_header_with_umi(header, umi):
    """This function inserts a UMI after the '@' symbol, making the downstream analysis easier
"""
    col=header.split(" ")
    # return "@" + umi + "_" + header[1:] + " " + head[2]
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
    parser = argparse.ArgumentParser(description='A script to identify reads in a UMI small RNA-seq fastq file',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input',
                        help='the input fastq file', required=True)
    parser.add_argument('-o', '--output',
                        help='the output fastq file', required=True)
    parser.add_argument('-d', '--pcr-duplicate',
                        help='the output fastq file containing PCR duplicates', required=True)
    parser.add_argument('-v', '--verbose', help='Also include detailed run info', action="store_true")
    parser.add_argument('-5', '--umi-pattern-5',
                        help='Set the UMI pattern at the 5\' end',
                        default='NNNCGANNNTACNNN')
    parser.add_argument('-3', '--umi-pattern-3',
                        help='Set the UMI pattern at the 3\' end',
                        default='NNNGTCNNNTAGNNN')
    parser.add_argument('-D', '--debug',
                        help='Turn on debugging mode', action="store_true")
    parser.add_argument('-q', '--quality',
                        help='Quality (phred quality score) cutoff for UMI. Default is 13, \
                        that is UMI with qualities >= 13 will be kept. This program assumes \
                        the phred quality scores in the fastq file are using sanger format',
                        required=False, type=int, default=13)
    # Sanger format can encode a Phred quality score from 0 to 93 using ASCII 33 to 126 (although in raw read data
    # the Phred quality score rarely exceeds 60, higher scores are possible in assemblies or read maps). Also used
    # in SAM format.[4] Coming to the end of February 2011, Illumina's newest version (1.8) of their pipeline CASAVA
    # will directly produce fastq in Sanger format, according to the announcement on seqanswers.com forum.[5]
    # These variables do not change during one run. They are assgined values once in the main() function
    global DEBUG
    global qc
    # global umi_len
    global c                              # Read count
    global umi_pat_5
    global umi_pat_3
    c = 0
    args = parser.parse_args()
    DEBUG = args.debug
    if DEBUG:
        print >>sys.stderr, "Debugging mode is on"
    # Quality cutoff
    qc = chr(args.quality + 33)
    print >>sys.stderr, "Quality cutoff in ASCII: " + qc
    fn = args.input
    verbose = args.verbose
    umi_pat_5 = args.umi_pattern_5
    umi_pat_3 = args.umi_pattern_3
    
    out = open(args.output, "w")

    n_drop_due_to_low_quality = 0
    # n_proper = 0

    # This stores the stats for each read
    stats = SraRunStats()
    if is_gzipped(fn):
        f = gzip.open(fn)
        print >>sys.stderr, "Input is gzipped."
    else:
        f = open(fn)
    while True:
        c += 1
        if c%4==1:
            r_name = f.readline().strip()
            if not r_name:
                break
        elif c % 4 == 2:
            r_seq = f.readline().strip()
        elif c % 4 == 3:
            r_info = f.readline().strip()
        else:
            r_qual = f.readline().strip()
            r = Read(r_name, r_seq, r_info, r_qual, "r1")
            r_name_proc, r_seq_proc, r_info_proc, r_qual_proc, r_bc = process_sra_read(r, stats)
            # print "#" + r1_name_proc + "#" + str(r1_name_proc == "") + "#" + str(r2_name_proc=="")
            # print "#" + r2_name_proc
            if r_name_proc == "":
                n_drop_due_to_low_quality += 1
                # print "Dropped read: ???"
            else:
                r_name_proc = get_header_with_umi(r_name_proc, r_bc)
                stats["n_with_proper_umi"] += 1
                print >>out, r_name_proc
                print >>out, r_seq_proc
                print >>out, r_info_proc
                print >>out, r_qual_proc

    f.close()
    print >>sys.stderr, ""
    print >>sys.stderr, "Stats: "
    # stats.summary()
        
    print >>sys.stderr, ""
    print >>sys.stderr, "Additional reads dropped because of low quality:\t" + str(n_drop_due_to_low_quality)
    print >>sys.stderr, "Final proper read pairs:\t" + str(stats["n_with_proper_umi"])
        
    print2("")
    if verbose:

        print2("=" * 80)
        print2(i + " padding usage")
        stats[i].padding_usage()
            
        print2("=" * 80)
        print2(i + " ligation bias")
        stats[i].ligation_bias()

if __name__ == "__main__":
    main()
