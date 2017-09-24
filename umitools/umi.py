#!/usr/bin/env python3

## TODO: put the classes for RNA-seq here as well


class SraUmiInfo:
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

