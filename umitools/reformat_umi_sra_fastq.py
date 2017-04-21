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
    def __init__(this, umi_pat5, umi_pat3):
        # Note that both should be lists so that it allows the two or more 5' UMI patterns
        this.pat5 = umi_pat5
        this.pat3 = umi_pat3

    def __str__(this):
        t = []
        t.append("5' UMI pattern: %s" % (",".join(this.pat5)))
        t.append("3' UMI pattern: %s" % (",".join(this.pat3)))
        return "\n".join(t)

    def check_umi_pat5_one(this, pat, r):
        '''it checks the reads against one 5' UMI pattern and returns
        the number of mismatches in fixed portion of the UMI and the
        nt in the variable portion of the read
'''
        bc = []
        mm = 0
        assert len(r) > len(pat), "Read is shorter than 5' UMI pattern: %s" % (r, )
        for i in range(len(pat)):
            if pat[i] == "N":
                bc.append(r[i])
            else:
                if r[i] != pat[i]:
                    mm += 1
        return [mm, "".join(bc)]

    def check_umi_pat3_one(this, pat, r):
        '''it checks the reads against one 3' UMI pattern and returns
        the number of mismatches in fixed portion of the UMI and the
        nt in the variable portion of the read
'''
        bc = []
        mm = 0
        offset = len(r) - len(pat)
        assert len(r) > len(pat), "Read is shorter than 3' UMI pattern: %s" % (r, )
        for i in range(len(pat)):
            if pat[i] == "N":
                bc.append(r[i+offset])
            else:
                if r[i+offset] != pat[i]:
                    mm += 1
        return [mm, "".join(bc)]

    def check(this, r):
        '''This function returns the total number of mismatches at 5' and
        3' end UMIs, together with the barcode (the variable portion of
        the UMI)
        '''
        mm5 = 100
        bc5 = ""
        for p in this.pat5:
            tmp_mm5, tmp_bc5 = this.check_umi_pat5_one(p, r)
            if tmp_mm5 < mm5:
                mm5 = tmp_mm5
                bc5 = tmp_bc5

        mm3 = 100
        bc3 = ""
        for p in this.pat3:
            tmp_mm3, tmp_bc3 = this.check_umi_pat3_one(p, r)
            if tmp_mm3 < mm3:
                mm3 = tmp_mm3
                bc3 = tmp_bc3

        # print "5' UMI mismatch %d, barcode %s" % (mm5, bc5)
        # print "3' UMI mismatch %d, barcode %s" % (mm3, bc3)
        return [mm5 + mm3, bc5 + bc3]
        
    def extract_insert(this, r):
        '''This function returns the insert part of reads. In other words, this 
        function extracts the actual small RNAs. It assumes all possible 5' UMI
        has the same lengths and all possible 3' UMIs have the same lengths
        '''
        return r[len(this.pat5[0]): len(r) - len(this.pat3[0])]

    def extract_insert_qual(this, r_qual):
        '''This function returns the qualities of the small RNA bases. It should be 
        used together with extrat_insert()
        '''
        return r_qual[len(this.pat5[0]): len(r_qual) - len(this.pat3[0])]


class SraRunStats():
    """It records the info for run (for small RNA-seq data)
"""
    def __init__(self):
        self.stats = {}
        self.stats["n_with_proper_umi"] = 0
        self.stats["n_non_duplicate"] = 0
        self.stats["n_duplpicate"] = 0        
        self.stats["n_without_proper_umi"] = 0
        self.stats["n_read"] = 0        

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
                        help='Set the UMI pattern at the 5\' end. Use ACGT for fixed nt and N for variable nt in UMI. If there are multiple patterns, separate them using comma',
                        default='NNNCGANNNTACNNN,NNNATCNNNAGTNNN')
    parser.add_argument('-3', '--umi-pattern-3',
                        help='Set the UMI pattern at the 3\' end. Use ACGT for fixed nt and N for variable nt in UMI. If there are multiple patterns, separate them using comma',
                        default='NNNGTCNNNTAGNNN')
    # Quality cutoff is not necessary as the tools for trimming 3' small RNA-seq adapters
    # already do this
    # 
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
