#!/usr/bin/env python3

import sys
from struct import unpack
import operator

__author__ = "Yu Fu"
__license__ = "GPLv3"


def print2(a):
    # print >>sys.stderr, a
    sys.stderr.write(str(a) + "\n")




class RsqRunStats():
    """A class for stats used in RNA-seq data
"""
    def __init__(self):
        self.stats = {}
        self.stats["n_read"] = 0
        self.stats["all_umi_locator"] = {}
        self.stats["n_without_locator"] = 0
        self.stats["n_with_locator"] = 0
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
        print2("Total:\t" + str(self["n_read"]))
        print2("Reads w/o locator:\t" + str(self["n_without_locator"]))
        print2("Reads w/ locator:\t" + str(self["n_with_locator"]))
        print2("-" * 80)
        print2("Reads w/ N's in UMI:\t" + str(self["n_with_ambiguous_umi"]))
        print2("Reads w/ wrong padding nt:\t" + 
               str(self["n_with_wrong_padding"]))
        print2("Reads w/ low-quality UMI:\t" + str(self["n_bad_quality_umi"]))
        print2("-" * 80)
        print2("Reads w/ proper UMI:\t" + str(self["n_good_reads"]))
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


class SraUmiInfo:
    '''A class for UMI information in small RNA-seq'''
    def __init__(this, umi_pat5, umi_pat3):
        '''Note that both should be lists so that it allows
        the two or more 5' UMI patterns'''
        this.pat5 = umi_pat5
        this.pat3 = umi_pat3

    def __str__(this):
        t = []
        t.append("5' UMI pattern: %s" % (",".join(this.pat5)))
        t.append("3' UMI pattern: %s" % (",".join(this.pat3)))
        return "\n".join(t)

    def check_umi_pat5_one(this, pat, r):
        '''it checks one read against one 5' UMI pattern and returns
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
        '''it checks one read against one 3' UMI pattern and returns
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
        '''This function returns the total number of mismatches of UMI locators at 5' and
        3' end of UMIs, together with the barcode (the variable portion of
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

#     def umi_locator_check(this, r):
#         '''Checks how many errors there are in the UMI locator portion of the read
# '''
#         bc = []
#         mm = 0
#         for i in range(len(pat)):
#             if pat[i] != "N":
#                 if r[i] != pat[i]:
#                     mm += 1
#         return [mm, "".join(bc)]

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
        self.stats["n_duplicate"] = 0
        self.stats["n_without_proper_umi"] = 0
        self.stats["n_read"] = 0

    def __getitem__(self, key):
        return self.stats[key]
    
    def __setitem__(self, key, value):
        self.stats[key] = value

    def report(self):
        sys.stderr.write("\n")
        sys.stderr.write("Stats: \n")
        sys.stderr.write("Total input reads:\t" + str(self["n_read"]) + "\n")
        sys.stderr.write("Reads dropped due to improper UMI:\t" +
                         str(self["n_without_proper_umi"]) + "\n")
        sys.stderr.write("Final proper read:\t" +
                         str(self["n_with_proper_umi"]) + "\n")
        sys.stderr.write("\tReads that are duplicates:\t" +
                         str(self["n_duplicate"]) + "\n")
        sys.stderr.write("\tReads that are non-duplicates:\t" +
                         str(self["n_non_duplicate"]) + "\n")


class SraRead():
    """ Informtion for one read including r_name, r_seq, 
    r_info, r_qual, mate and so on...
"""
    def __init__(self, r_name, r_seq, r_info, r_qual, mate):
        self.r_name = r_name
        self.r_seq = r_seq
        self.r_info = r_info
        self.r_qual = r_qual
        self.mate = mate

    def __str__(self):
        ret = "\n".join((self.r_name, self.r_seq, self.r_info, self.r_qual))
        ret = ret + "\n"
        return ret


def is_gzipped(filename):
    # 1F 8B 08 00 / gz magic number
    # 1F 8B is the magic number. The 08 and 00 are not always.
    # In python 2 magic = (b'\x1f', b'\x8b', b'\x08', b'\x00') is fine
    # In python 3, I need to use magic = (b'\x1f', b'\x8b', b'\x08', b'\x00')
    # magic = (b'\x1f', b'\x8b', b'\x08', b'\x00')
    magic = (b'\x1f', b'\x8b')
    with open(filename, 'rb') as handle:
        s = unpack('cc', handle.read(2))
        return s == magic


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


class bcolors:
    '''From https://stackoverflow.com/questions/287871/
    print-in-terminal-with-colors-using-python
    '''
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

