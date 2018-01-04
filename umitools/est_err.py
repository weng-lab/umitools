#!/usr/bin/env python3

# Estimate the error rates of short reads using the fixed portion of reads
# in UMI small RNA-seq

# TODO: Do I need to do this for UMI RNA-seq?

import umitools.umi
import argparse
import gzip


def check_one_umi_locator(pat, r):
    pos_matched = []
    for i in range(len(pat)):
        if pat[i] != "N":
            if pat[i] != r[i]:
                pos_matched.append(False)
            else:
                pos_matched.append(True)

    return(pos_matched)


def check_umi_locator(pats, r):
    '''This function iterates all UMI patterns found in pats and
    report positions of matches for the best UMI. The return value
    is a list of boolean values indicating if a position is a match
'''
    max_matches = 0
    pos_matched = []
    for pat in pats:
        tmp = check_one_umi_locator(pat, r)
        if(sum(tmp) > max_matches):
            max_matches = max_matches
            pos_matched = tmp
    return(pos_matched)

def test():
    r = 'GAAAACCAGAGTCATTCTCTGGCAGATCAACTGTTCTCTGCATTTGATGTCAAGTAGCGT'
    umi5_len = len(umi_pat5[0])
    r5 = r[0:umi5_len]
    pos_mat5 = check_umi_locator(umi_pat5, r5)

    umi3_len = len(umi_pat3[0])
    r3 = r[len(r)-umi3_len:]
    pos_mat3 = check_umi_locator(umi_pat3, r3)

    pos_mat = pos_mat5 + pos_mat3
    print("Total number of matches in UMI locator: {}".
          format(sum(pos_mat)))

    print("Total number of mismatches in UMI locator: {}".
          format(len(pos_mat) - sum(pos_mat)))

    # One mismatch in the locator
    r = "NNNCGANNNTAGNNN"
    print(check_umi_locator(umi_pat5, r))

    # Two mismatches in the locator
    r = "NNNGGANNNTAGNNN"
    print(check_umi_locator(umi_pat5, r))

    # For UMI locator at the 3' end

def main():

    desc = "Is A script that estimates the sequencing error rates"
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', help="Input file name")
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

    args = parser.parse_args()
    umi_pat5 = args.umi_pattern_5.split(",")
    umi_pat3 = args.umi_pattern_3.split(",")
    if umi.is_gzipped(args.f):
        print("Input is gzipped")
        f = gzip.open(args.f, "rt")
    else:
        f = open(args.f)

    M = []                      # How many reads are matches at each position
    tmp_l = len(umi_pat5[0]) - umi_pat5[0].count("N") + \
            len(umi_pat3[0]) - umi_pat3[0].count("N")
    for i in range(tmp_l):
        M.append(0)

    # test()
    c = 0
    total_read = 0
    total_valid_read = 0
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
            # r = umi.SraRead(r_name, r_seq, r_info, r_qual, "r1")
            # r_name_proc, r_seq_proc, r_info_proc, \
            #     r_qual_proc, r_bc = process_sra_read(r, stats, ui)
            
            r = r_seq
            umi5_len = len(umi_pat5[0])
            r5 = r[0:umi5_len]
            pos_mat5 = check_umi_locator(umi_pat5, r5)

            umi3_len = len(umi_pat3[0])
            r3 = r[len(r)-umi3_len:]
            pos_mat3 = check_umi_locator(umi_pat3, r3)

            pos_mat = pos_mat5 + pos_mat3
            total_read += 1
            # print("Total number of matches in UMI locator: {}".
            #       format(sum(pos_mat)))
            # print("Total number of mismatches in UMI locator: {}".
            #       format(len(pos_mat) - sum(pos_mat)))

            if sum(pos_mat) >= len(pos_mat) - 1:
                # The valid reads should be reported instead of total reads
                # because some reads are just garbage
                total_valid_read += 1
                for i in range(len(pos_mat)):
                    if pos_mat[i]:
                        M[i] += 1

    print("UMI_locator_pos\tmatch\ttotal")
    for i in range(len(M)):
        print("{}\t{}\t{}".format(i, M[i], total_valid_read))

if __name__ == "__main__":
    main()
    
