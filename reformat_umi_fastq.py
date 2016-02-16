#!/usr/bin/env python

import gzip
import sys
import argparse

def process_read(r_name, r_seq, r_info, r_qual, mate, stats):
    """This function processes one read. mate can be "r1" or "r2". Returns 4 empty strings if the read is dropped
"""
    ret_name = ""
    ret_seq = ""
    ret_info = ""
    ret_qual = ""
    ret_bc = ""
    if DEBUG == True:
        print '-' * 80
        print mate
        print r_name
        print "Original qual:\t" + r_qual
        print "Original read:\t" + r_seq
        print "Locator:\t" + ' ' * 5 + r_seq[umi_len : umi_len + umi_locator_len]
        print "UMI:\t\t" + r_seq[0: umi_len]
        print "What's left:\t" + ' ' * 9 + r_seq[umi_len + umi_locator_len + umi_downstream_len: ]
    if r_seq[umi_len : umi_len + umi_locator_len] == umi_locator:
        stats[mate]["n_with_locator"] += 1
        my_umi = r_seq[0:umi_len]
        my_umi_qual = r_qual[0:umi_len]
        if my_umi.find("N") == -1:
            if r_seq[umi_len + umi_locator_len] == umi_downstream:
                if is_good_phred(my_umi_qual, qc):
                    my_seq = r_seq[umi_len + umi_locator_len + umi_downstream_len: ]
                    stats[mate]["n_good_reads"] += 1
                    # I do not modify read name here,
                    # since we need to store barcodes from both reads and put the
                    # concatenated barcode into both reads to make downstream analysis easier
                    # ret_name = get_header_with_umi(r_name, my_umi)
                    ret_name = r_name
                    ret_seq =  my_seq
                    ret_info = r_info
                    ret_qual = r_qual[umi_len + umi_locator_len + umi_downstream_len: ]
                    ret_bc = my_umi
                else:
                    stats[mate]["n_bad_quality_umi"] += 1
                    if DEBUG:
                        print "*Found one with low quality UMI:\t" + r_seq
            else:
                stats[mate]["n_with_wrong_padding"] += 1
                if DEBUG:
                    print "*Found one with wrong padding nucleotide (A/C/G):\t" + r_seq
        else:
            if DEBUG == True:
                print "*Found one with N's in UMI:\t" + r_seq
            stats[mate]["n_with_ambiguous_umi"] += 1
    else:
        if DEBUG == True:
            print "*Found one without locator:\t" + r_seq
        stats[mate]["n_without_locator"] += 1
    return ret_name, ret_seq, ret_info, ret_qual, ret_bc

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

# print phred_checker("BBC", "B")
# print phred_checker("ABC", "B")

parser = argparse.ArgumentParser(description='A script to reformat r1 reads in a UMI fastq file so that the name of each record contains the UMI')
parser.add_argument('-l', '--left', help='the input fastq file for r1. If you want to pipe in the input, you can use /dev/stdin', required=True)
parser.add_argument('-r', '--right', help='the input fastq file for r2. If you want to pipe in the input, you can use /dev/stdin', required=True)
parser.add_argument('-L', '--left-out', help='the output fastq file for r1', required=True)
parser.add_argument('-R', '--right-out', help='the output fastq file for r2', required=True)
parser.add_argument('-q', '--quality', help='Quality (phred quality score) cutoff for UMI. Default is 20, that is UMI with qualities >= 20 will be kept. This program assumes the phred quality scores in the fastq file are using sanger format', required=False, type=int, default=20)
parser.add_argument('-D', '--debug', help='Turn on debugging mode', action="store_true")

# Sanger format can encode a Phred quality score from 0 to 93 using ASCII 33 to 126 (although in raw read data the Phred quality score rarely exceeds 60, higher scores are possible in assemblies or read maps). Also used in SAM format.[4] Coming to the end of February 2011, Illumina's newest version (1.8) of their pipeline CASAVA will directly produce fastq in Sanger format, according to the announcement on seqanswers.com forum.[5]

args = parser.parse_args()
DEBUG=args.debug
if DEBUG:
    print >>sys.stderr, "Debugging mode is on"
# Quality cutoff
qc = chr(args.quality + 33)
print >>sys.stderr, "Quality cutoff in ASCII:\t" + qc
fn1 = args.left
fn2 = args.right

out1 = open(args.left_out, "w")
out2 = open(args.right_out, "w")

c = 0

umi_len = 5
umi_locator = 'GGG'
umi_locator_len = len(umi_locator)

# Trim one nucleotide after the GGG
umi_downstream = 'T'
umi_downstream_len = len(umi_downstream)

n_additional_drop_due_to_mate = 0
n_proper_pair = 0

# This stores the stats for each read
stats = {}
for i in ("r1", "r2"):
    # Those reads without or with GGG
    stats[i] = {}
    stats[i]["n_without_locator"] = 0
    stats[i]["n_with_locator"] = 0
    # Those reads with N's before GGG
    stats[i]["n_with_ambiguous_umi"] = 0
    # Those with A/C/G after GGG
    stats[i]["n_with_wrong_padding"] = 0
    # Those having bad quality UMI w/ GGG and T and w/o N's
    stats[i]["n_bad_quality_umi"] = 0
    # The rest
    stats[i]["n_good_reads"] = 0
f1 = open(fn1)
f2 = open(fn2)

while True:
    c += 1
    if c%4==1:
        r1_name = f1.readline().strip()
        r2_name = f2.readline().strip()
        if not r1_name:
            break
    elif c%4==2:
        r1_seq = f1.readline().strip()
        r2_seq = f2.readline().strip()
    elif c%4==3:
        r1_info = f1.readline().strip()
        r2_info = f2.readline().strip()
    else:
        r1_qual = f1.readline().strip()
        r2_qual = f2.readline().strip()
        r1_name_proc, r1_seq_proc, r1_info_proc, r1_qual_proc, r1_bc = process_read(r1_name, r1_seq, r1_info, r1_qual, "r1", stats)
        r2_name_proc, r2_seq_proc, r2_info_proc, r2_qual_proc, r2_bc = process_read(r2_name, r2_seq, r2_info, r2_qual, "r2", stats)
        # print "#" + r1_name_proc + "#" + str(r1_name_proc == "") + "#" + str(r2_name_proc=="")
        # print "#" + r2_name_proc
        if (r1_name_proc == "" and r2_name_proc != "") or (r2_name_proc == "" and r1_name_proc != ""):
            n_additional_drop_due_to_mate += 1
        elif r1_name_proc != "" and r2_name_proc != "":
            r1_name_proc = get_header_with_umi(r1_name_proc, r1_bc + r2_bc)
            r2_name_proc = get_header_with_umi(r2_name_proc, r1_bc + r2_bc)            
            n_proper_pair += 1
            print >>out1, r1_name_proc
            print >>out1, r1_seq_proc
            print >>out1, r1_info_proc
            print >>out1, r1_qual_proc
            print >>out2, r2_name_proc
            print >>out2, r2_seq_proc
            print >>out2, r2_info_proc
            print >>out2, r2_qual_proc

f1.close()
f2.close()
# out1.close()
# out2.close()
for i in ("r1", "r2"):
    print >>sys.stderr, "-" * 80
    print >>sys.stderr, i
    print >>sys.stderr, "Total:\t" + str((c-1)/4)
    print >>sys.stderr, "Reads w/o locator:\t" + str(stats[i]["n_without_locator"])
    print >>sys.stderr, "Reads w/ locator:\t" + str(stats[i]["n_with_locator"])
    print >>sys.stderr, "-" * 80
    print >>sys.stderr, "Reads w/ N's in UMI:\t" + str(stats[i]["n_with_ambiguous_umi"])
    print >>sys.stderr, "Reads w/ wrong padding nt (A/C/G):\t" + str(stats[i]["n_with_wrong_padding"])
    print >>sys.stderr, "Reads w/ low-quality UMI:\t" + str(stats[i]["n_bad_quality_umi"])
    print >>sys.stderr, "-" * 80
    print >>sys.stderr, "Reads w/ proper UMI:\t" + str(stats[i]["n_good_reads"])
    print >>sys.stderr, "-" * 80
print >>sys.stderr, ""
print >>sys.stderr, "Additional reads dropped because its mate is dropped:\t" + str(n_additional_drop_due_to_mate)
print >>sys.stderr, "Final proper read pairs:\t" + str(n_proper_pair)
