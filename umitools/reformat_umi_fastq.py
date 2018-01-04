#!/usr/bin/env python3

import gzip
import sys
import argparse
import re
from struct import unpack
import umitools.umi as umi

__author__ = "Yu Fu"
__license__ = "GPLv3"


class RsqRead():
    """Informtion for one RNA-seq read including r_name, r_seq, r_info, r_qual, mate
    and so on...
    """
    def __init__(self, r_name, r_seq, r_info, r_qual, mate):
        self.r_name = r_name
        self.r_seq = r_seq
        self.r_info = r_info
        self.r_qual = r_qual
        self.mate = mate
        self.my_umi_locator = r_seq[umi_len : umi_len + umi_locator_len]
        self.my_umi = my_umi = r_seq[0:umi_len]
        self.my_umi_qual = r_qual[0:umi_len]


def process_read(read, stats):
    """This function processes one read. The "mate" option can be "r1" or "r2".
    Returns 4 empty strings if the read is dropped
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
        print('-' * 80)
        print(mate)
        print(r_name)
        print("Original qual:\t" + r_qual)
        print("Original read:\t" + r_seq)
        print("Locator:\t" + ' ' * 5 + 
              r_seq[umi_len: umi_len + umi_locator_len])
        print("UMI:\t\t" + r_seq[0: umi_len])
        print("What's left:\t" + ' ' * 9 + 
              r_seq[umi_len + umi_locator_len + 
                    umi_downstream_len:])
    my_umi_locator = read.my_umi_locator
    my_umi = read.my_umi
    if my_umi_locator in stats[mate]["all_umi_locator"]:
        stats[mate]["all_umi_locator"][my_umi_locator] += 1
    else:
        stats[mate]["all_umi_locator"][my_umi_locator] = 1
    if my_umi_locator in umi_locators:
        stats[mate]["n_with_locator"] += 1
        my_umi_qual = read.my_umi_qual  # r_qual[0:umi_len]
        if my_umi.find("N") == -1:
            padding = r_seq[umi_len + umi_locator_len]
            if padding in umi_downstream:
                stats[mate]["padding"][padding] += 1
                if umi.is_good_phred(my_umi_qual, qc):
                    my_seq = r_seq[umi_len + 
                                   umi_locator_len + umi_downstream_len:]
                    stats[mate]["n_good_reads"] += 1
                    stats[mate]["ligation"][my_seq[0]] += 1
                    # I do not modify read name here,
                    # since we need to store barcodes from both reads and put
                    # the concatenated barcode into both reads to make
                    # downstream analysis easier
                    # ret_name = get_header_with_umi(r_name, my_umi)
                    ret_name = r_name
                    ret_seq = my_seq
                    ret_info = r_info
                    ret_qual = r_qual[umi_len + umi_locator_len + 
                                      umi_downstream_len:]
                    ret_bc = my_umi
                    if my_umi_locator in stats[mate]["good_umi_locator"]:
                        stats[mate]["good_umi_locator"][my_umi_locator] += 1
                    else:
                        stats[mate]["good_umi_locator"][my_umi_locator] = 1
                else:
                    stats[mate]["n_bad_quality_umi"] += 1
                    if DEBUG:
                        print("*Found one with low quality UMI:\t" + r_seq)
            else:
                stats[mate]["n_with_wrong_padding"] += 1
                if DEBUG:
                    print("*Found one with wrong padding nucleotide (A/C/G):\t" 
                          + r_seq)
        else:
            if DEBUG:
                print("*Found one with N's in UMI:\t" + r_seq)
            stats[mate]["n_with_ambiguous_umi"] += 1
    else:
        if DEBUG:
            print("*Found one without locator:\t" + r_seq)
        stats[mate]["n_without_locator"] += 1
    return ret_name, ret_seq, ret_info, ret_qual, ret_bc


def get_header_with_umi(header, umi):
    '''This function inserts a UMI after the '@' symbol, making the
    downstream analysis easier
    '''
    col = header.split(" ")
    # return "@" + umi + "_" + header[1:] + " " + head[2]
    return col[0] + "_" + umi + " " + col[1]


def main():
    parser = argparse.ArgumentParser(description='''
A script to reformat
reads in a UMI fastq file so that the name of each record contains the 
UMI. This script is also known as umitools extract.''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-l', '--left', help='the input fastq file for r1.',
                        required=True)
    parser.add_argument('-r', '--right', help='the input fastq file for r2.',
                        required=True)
    parser.add_argument('-L', '--left-out',
                        help='the output fastq file for r1',
                        required=True)
    parser.add_argument('-R', '--right-out',
                        help='the output fastq file for r2',
                        required=True)
    parser.add_argument('-v', '--verbose', help='Also include detailed stats for \
    UMI and padding usage', action="store_true")
    parser.add_argument('--umi-locator',
                        help="Set the UMI locators. If you have multiple, "
                        "separate them by comma. e.g. GGG,TCA,ATC",
                        default='GGG,TCA,ATC')
    parser.add_argument('--umi-padding',
                        help='''
Set the nucleotide (for preventing ligation bias)
after the UMI locators. If you have multiple, separate
them by comma. e.g. A,C,G,T. The quality for this nt is
sometimes low, so the default is all possible four
nucleotides''', default='A,C,G,T,N')
    # An alternative way to specify the UMI pattern
    parser.add_argument('--umi-pattern', help='Set the UMI patterns.')

    parser.add_argument('-q', '--quality',
                        help='Quality (phred quality score) cutoff for UMI.'
                        'Default is 13, that is UMI with qualities >= 13 will'
                        'be kept. This program assumes the phred quality scores'
                        'in the fastq file are using sanger format',
                        required=False, type=int, default=13)

    parser.add_argument('-D', '--debug', help='Turn on debugging mode',
                        action="store_true")

    # For details on Sanger format for phred scores, check out
    # https://en.wikipedia.org/wiki/FASTQ_format
    global DEBUG
    global qc
    global umi_len
    global umi_locators
    global umi_locator_len
    global umi_downstream
    global umi_downstream_len

    c = 0
    args = parser.parse_args()
    DEBUG = args.debug
    if DEBUG:
        umi.print2("Debugging mode is on")
    # Quality cutoff
    qc = chr(args.quality + 33)
    umi.print2("Quality cutoff in ASCII: " + qc)
    fn1 = args.left
    fn2 = args.right
    verbose = args.verbose
    if re.search("\.gz|\.gzip", args.left_out):
        out1 = gzip.open(args.left_out, "wt", compresslevel=4)
    else:
        out1 = open(args.left_out, "w")
    if re.search("\.gz|\.gzip", args.right_out):
        out2 = gzip.open(args.right_out, "wt", compresslevel=4)
    else:
        out2 = open(args.right_out, "w")

    umi_len = 5
    umi_locators = {}
    tmp = args.umi_locator.split(",")
    for i in tmp:
        umi_locators[i] = True
    # umi_locators = ['GGG', 'GAT', 'TGA']
    # umi.print2 >>sys.stderr, "UMI locators: ",
    sys.stderr.write("UMI locators: ",)
    for i in umi_locators:
        sys.stderr.write(i)
    # print >>sys.stderr, ""
    sys.stderr.write("\n")
    umi_locator_len = len(tmp[0])
    # Trim one nucleotide after the GGG
    # umi_downstream = 'T'
    umi_downstream = {}
    tmp = args.umi_padding.split(",")
    for i in tmp:
        umi_downstream[i] = True
    umi_downstream_len = len(tmp[0])
    sys.stderr.write("UMI padding: ")
    for i in umi_downstream:
        sys.stderr.write(i)
    sys.stderr.write("\n")

    n_additional_drop_due_to_mate = 0
    n_proper_pair = 0

    # This stores the stats for each read
    stats = {}
    for i in ("r1", "r2"):
        stats[i] = umi.RsqRunStats()
    f1 = open(fn1)
    f2 = open(fn2)
    if umi.is_gzipped(fn1):
        f1 = gzip.open(fn1, "rt")
        umi.print2("r1 input is gzipped.")
    else:
        f1 = open(fn1)
    if umi.is_gzipped(fn2):
        f2 = gzip.open(fn2, "rt")
        umi.print2("r2 input is gzipped.")
    else:
        f2 = open(fn2)
    while True:
        c += 1
        if c % 4 == 1:
            r1_name = f1.readline().strip()
            r2_name = f2.readline().strip()
            if not r1_name:
                break
        elif c % 4 == 2:
            r1_seq = f1.readline().strip()
            r2_seq = f2.readline().strip()
            stats["r1"]["n_read"] += 1
            stats["r2"]["n_read"] += 1
        elif c % 4 == 3:
            r1_info = f1.readline().strip()
            r2_info = f2.readline().strip()
        else:
            r1_qual = f1.readline().strip()
            r2_qual = f2.readline().strip()
            r1 = RsqRead(r1_name, r1_seq, r1_info, r1_qual, "r1")
            r2 = RsqRead(r2_name, r2_seq, r2_info, r2_qual, "r2")
            r1_name_proc, r1_seq_proc, r1_info_proc, r1_qual_proc, r1_bc = process_read(r1, stats)
            r2_name_proc, r2_seq_proc, r2_info_proc, r2_qual_proc, r2_bc = process_read(r2, stats)
            # print "#" + r1_name_proc + "#" + str(r1_name_proc == "") + "#" + str(r2_name_proc=="")
            # print "#" + r2_name_proc
            if (r1_name_proc == "" and r2_name_proc != "") or (r2_name_proc == "" and r1_name_proc != ""):
                n_additional_drop_due_to_mate += 1
            elif r1_name_proc != "" and r2_name_proc != "":
                r1_name_proc = get_header_with_umi(r1_name_proc, r1_bc + r2_bc)
                r2_name_proc = get_header_with_umi(r2_name_proc, r1_bc + r2_bc)            
                n_proper_pair += 1
                print(r1_name_proc, file=out1)
                print(r1_seq_proc, file=out1)
                print(r1_info_proc, file=out1)
                print(r1_qual_proc, file=out1)
                print(r2_name_proc, file=out2)
                print(r2_seq_proc, file=out2)
                print(r2_info_proc, file=out2)
                print(r2_qual_proc, file=out2)

    f1.close()
    f2.close()
    # out1.close()
    # out2.close()
    # print >>sys.stderr, ""
    sys.stderr.write("\n")
    for i in ("r1", "r2"):
        # print >>sys.stderr, "Stats for " + i
        sys.stderr.write("Stats for " + i + "\n")
        stats[i].summary()
        
    # print >>sys.stderr, ""
    umi.print2("")

    # print >>sys.stderr, "Additional reads dropped because its mate is dropped:\t" + str(n_additional_drop_due_to_mate)
    umi.print2("Additional reads dropped because its mate is dropped:\t" + 
               str(n_additional_drop_due_to_mate))
    # print >>sys.stderr, "Final proper read pairs:\t" + str(n_proper_pair)
    umi.print2("Final proper read pairs:\t" + str(n_proper_pair))
    
    umi.print2("")
    if verbose:
        for i in ("r1", "r2"):
            umi.print2("=" * 80)
            umi.print2(i + " UMI padding usage")
            stats[i].umi_padding_usage()

            umi.print2("=" * 80)
            umi.print2(i + " padding usage")
            stats[i].padding_usage()
            
            umi.print2("=" * 80)
            umi.print2(i + " ligation bias")
            stats[i].ligation_bias()


if __name__ == "__main__":
    main()
