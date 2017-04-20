#!/usr/bin/env python

import gzip
import sys
import argparse
from struct import unpack

__author__ = "Yu Fu"
__license__ = "GPLv3"


def print2(a):
    print >>sys.stderr, a

class UMIInfo:
    def transform_umi_to_pos(this, umi_pat):
        '''Transform the UMI pattern into the positions of N's and positions of 
non-N's, which makes the UMI pattern searching faster
'''
        a = zip(range(0, len(umi_pat)), list(umi_pat))
        ret_pos_n = [i[0] for i in a if i[1] == "N"]
        ret_pos_fixed = [i[0] for i in a if i[1] != "N"]
        ret_base_fixed = [i[1] for i in a if i[1] != "N"]
        return (ret_pos_n, ret_pos_fixed, ret_base_fixed)

    def __init__(this, umi_pat5, umi_pat3):
        this.pat5 = umi_pat5
        this.pat3 = umi_pat3
        this.pos_n5, this.pos_fixed5, this.base_fixed5 = this.transform_umi_to_pos(umi_pat5)
        this.pos_n3, this.pos_fixed3, this.base_fixed3 = this.transform_umi_to_pos(umi_pat3)            

    def __str__(this):
        t = []
        t.append("5' UMI pattern: %s" % (this.pat5))
        t.append("Positions of variable nt:")
        t.append(",".join([str(i) for i in this.pos_n5]))
        t.append("Postions of fixed nt:")
        t.append(",".join([str(i) for i in this.pos_fixed5]))
        t.append("Fixed nt:")
        t.append(",".join([str(i) for i in this.base_fixed5]))
        t.append("")
        t.append("3' UMI pattern: %s" % (this.pat3))
        t.append("Positions of variable nt:")
        t.append("this.pos_n3")
        t.append("Postions of fixed nt:")
        t.append(",".join([str(i) for i in this.pos_fixed3]))
        t.append("Fixed nt:")
        t.append(",".join([str(i) for i in this.base_fixed3]))
        return "\n".join(t)

    def n_mismatch(this, r):
        '''This function returns of the numbers of mismatches in the fixed portion of
5' UMI and 3' UMI'''
        a = [r[i] for i in this.pos_fixed5]
        b = this.base_fixed5
        assert len(a) == len(b)
        ret5 = 0
        for i in range(len(a)):
            if a[i] != b[i]:
                ret5 += 1

        ret3 = 0
        offset = len(r) - len(this.pat3)
        a = [r[i+offset] for i in this.pos_fixed3]
        b = this.base_fixed3
        assert len(a) == len(b)
        for i in range(len(a)):
            if a[i] != b[i]:
                ret3 += 1
        return (ret5, ret3)

    def extract_umi(this, r):
        '''This function extracts UMI info. In other words, this function extracts
        bases of reads that correspond to the N part of the UMI pattern
        '''
        offset = len(r) - len(this.pat3)
        p1 = "".join([r[i] for i in this.pos_n3])
        p2 = "".join([r[i+offset] for i in this.pos_n3])
        return p1 + p2

    def extract_insert(this, r):
        '''This function returns the insert part of reads. In other words, this 
        function extracts the actual small RNAs
        '''
        return r[len(this.pat5): len(r) - len(this.pat3)]

    def extract_insert_qual(this, r_qual):
        '''This function returns the qualities of the small RNA bases. It should be 
        used together with extrat_insert()
        '''
        return r_qual[len(this.pat5): len(r_qual) - len(this.pat3)]


class SraRunStats():
    """It records the info for run (for small RNA-seq data)
"""
    def __init__(self):
        self.stats = {}
        self.stats["n_with_proper_umi"] = 0
        self.stats["n_non_duplicate"] = 0
        self.stats["n_duplpicate"] = 0        
        self.stats["n_without_proper_umi"] = 0

    def __getitem__(self, key):
        return self.stats[key]
    def __setitem__(self, key, value):
        self.stats[key] = value


class Read():
    """ Informtion for one read including r_name, r_seq, 
    r_info, r_qual, mate and so on...
"""
    def __init__(self, r_name, r_seq, r_info, r_qual, mate):
        self.r_name = r_name
        self.r_seq = r_seq
        self.r_info = r_info
        self.r_qual = r_qual
        self.mate = mate


def is_gzipped(filename):
    # 1F 8B 08 00 / gz magic number
    magic = ('\x1f', '\x8b', '\x08', '\x00')

    with open(filename, 'rb') as handle:
        s = unpack('cccc', handle.read(4))
        return s == magic


# @profile
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
    mate = read.mate
    
    n_mm5, n_mm3 = ui.n_mismatch(r_seq)
    if n_mm5 + n_mm3 <= N_MISMATCH_ALLOWED_IN_UMI:
        ret_bc = ui.extract_umi(r_seq)
        ret_name = r_name
        ret_seq = ui.extract_insert(r_seq)
        ret_info = r_info
        ret_qual = ui.extract_insert_qual(r_qual)
        # print ret_seq
        # print ret_bc        
    return ret_name, ret_seq, ret_info, ret_qual, ret_bc

    
def get_header_with_umi(header, umi):
    """This function inserts a UMI after the '@' symbol, making the downstream 
    analysis easier
    """
    col=header.split(" ")
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
    desc = "A script to identify reads in a UMI small RNA-seq fastq file"
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input',
                        help='the input fastq file', required=True)
    parser.add_argument('-o', '--output',
                        help='the output fastq file', required=True)
    parser.add_argument('-d', '--pcr-duplicate',
                        help='the output fastq file containing PCR duplicates',
                        required=True) 
    parser.add_argument('--reads-with-improper-umi', default="",
                        help='the output fastq file containing reads with improper UMIs. \
                        The default is to throw away these reads. This is for debugging purposes',
                        required=False)
   
    parser.add_argument('-v', '--verbose',
                        help='Also include detailed run info',
                        action="store_true")
    parser.add_argument('-5', '--umi-pattern-5',
                        help='Set the UMI pattern at the 5\' end',
                        default='NNNCGANNNTACNNN')
    parser.add_argument('-3', '--umi-pattern-3',
                        help='Set the UMI pattern at the 3\' end',
                        default='NNNGTCNNNTAGNNN')
    # parser.add_argument('-D', '--debug',
    #                     help='Turn on debugging mode', action="store_true")
    # Quality cutoff is not necessary as the tools for trimming 3' small RNA-seq adapters
    # already do this
    # 
    # parser.add_argument('-q', '--quality',
    #                     help='Quality (phred quality score) cutoff for UMI. Default is 13, \
    #                     that is UMI with qualities >= 13 will be kept. This \
    #                     program assumes the phred quality scores in the fastq \
    #                     file are using sanger format',
    #                     required=False, type=int, default=13)
    # global DEBUG
    global qc
    global c # Read count
    global N_MISMATCH_ALLOWED_IN_UMI
    N_MISMATCH_ALLOWED_IN_UMI = 1
    c = 0
    args = parser.parse_args()
    # DEBUG = args.debug
    # if DEBUG:
    #     print >>sys.stderr, "Debugging mode is on"
    # Quality cutoff (This is not used for now as 3' adapter removal has already
    # taken care of this
    # qc = chr(args.quality + 33)
    # print >>sys.stderr, "Quality cutoff in ASCII: " + qc
    
    fn = args.input
    verbose = args.verbose
    umi_pat5 = args.umi_pattern_5
    umi_pat3 = args.umi_pattern_3

    ui = UMIInfo(umi_pat5, umi_pat3)

    print >>sys.stderr, "-" * 72
    print >>sys.stderr, "Summary of UMI patterns at 5' and 3' ends:"
    print >>sys.stderr, "-" * 72    
    print >>sys.stderr, ui
    print >>sys.stderr, "-" * 72
        
    out = open(args.output, "w")
    dup = open(args.pcr_duplicate, "w")
    if len(args.reads_with_improper_umi) != 0:
        imp = open(args.reads_with_improper_umi, "w")
        print args.reads_with_improper_umi
    n_drop_due_to_low_quality = 0
    # n_proper = 0

    stats = SraRunStats()
    if is_gzipped(fn):
        f = gzip.open(fn)
        print >>sys.stderr, "Input is gzipped."
    else:
        f = open(fn)

    ## Small RNA insert and umi should be unique; otherwise, it is a duplicate
    insert_umi = {}
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
            r_name_proc, r_seq_proc, r_info_proc, r_qual_proc, r_bc = process_sra_read(r, stats, ui)
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
    # stats.summary()
        
    print >>sys.stderr, ""
    print >>sys.stderr, "Reads dropped due to improper UMI:\t" + str(stats["n_without_proper_umi"])
    print >>sys.stderr, "Final proper read:\t" + str(stats["n_with_proper_umi"])
    print >>sys.stderr, "\tReads that are duplicates:\t" + str(stats["n_duplpicate"])
    print >>sys.stderr, "\tReads that are non-duplicates:\t" + str(stats["n_non_duplicate"])
        
    print2("")
    if verbose:
        pass

if __name__ == "__main__":
    main()
