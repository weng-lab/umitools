#!/usr/bin/env python3

import gzip
import sys
import argparse
from struct import unpack
import re

import umitools.umi as umi
import umitools.umi_graph as umi_graph

__author__ = "Yu Fu"
__license__ = "GPLv3"


def print2(a):
    # print >>sys.stderr, a
    sys.stderr.write(str(a) + "\n")

    
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

    if mm <= N_MISMATCH_ALLOWED_IN_UMI_LOCATOR:
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
    desc = '''A script to process reads in from UMI small RNA-seq. This script can handle
gzipped files transparently. This script is also known as
umitools extract_small.'''
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input',
                        help='the input fastq file.', required=True)
    parser.add_argument('-o', '--output',
                        help='''the output fastq file containing reads that
                        are not duplicates''', required=True)
    parser.add_argument('-d', '--pcr-duplicate',
                        help='The output fastq file containing PCR duplicates',
                        required=True) 
    parser.add_argument('--reads-with-improper-umi', default="",
                        help='The output fastq file containing reads with improper UMIs. \
                        The default is to throw away these reads. This is for \
                        debugging purposes',
                        required=False)
    parser.add_argument('-e', '--errors-allowed', type=int, default=0,
                        help="Setting it to >=1 allows errors in UMIs. \
                        Otherwise, no errors are allowed in UMIs.")
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

    # --force-graph is used for debugging. Ideally, when the network method
    # does not connect nodes at all, it should output the same results
    # with the 'unique' method.
    parser.add_argument('--force-graph', help=argparse.SUPPRESS,
                        action="store_true")

    parser.add_argument('--debug', help='More output for debugging', action="store_true")
    # Quality cutoff is not necessary as the tools for trimming 3' small RNA-seq adapters
    # already do this
    # parser.add_argument('-q', '--quality',
    #                     help='Quality (phred quality score) cutoff for UMI. Default is 13, \
    #                     that is UMI with qualities >= 13 will be kept. This \
    #                     program assumes the phred quality scores in the fastq \
    #                     file are using sanger format',
    #                     required=False, type=int, default=13)
    global N_MISMATCH_ALLOWED_IN_UMI_LOCATOR
    N_MISMATCH_ALLOWED_IN_UMI_LOCATOR = 1
    global DEBUG

    args = parser.parse_args()
    
    fn = args.input
    verbose = args.verbose
    umi_pat5 = args.umi_pattern_5.split(",")
    umi_pat3 = args.umi_pattern_3.split(",")
    umi_errors_allowed = args.errors_allowed
    DEBUG = args.debug
    FORCE_GRAPH = args.force_graph
    if FORCE_GRAPH:
        print2("Debugging: force the use of network method even " 
               "when no mismatch is allowed")
    ui = umi.SraUmiInfo(umi_pat5, umi_pat3)

    sys.stderr.write("-" * 72 + "\n")
    sys.stderr.write("Summary of UMI patterns at 5' and 3' ends:\n")
    sys.stderr.write("-" * 72 + "\n")
    sys.stderr.write(str(ui) + "\n")
    sys.stderr.write("-" * 72 + "\n")
    sys.stderr.write("Number of UMI errors allowed: {}".
                     format(umi_errors_allowed) + "\n")

    if re.search("\.gz|\.gzip", args.output):
        out = gzip.open(args.output, "wt", compresslevel=4)
    else:
        out = open(args.output, "w")
    if re.search("\.gz|\.gzip", args.pcr_duplicate):
        dup = gzip.open(args.pcr_duplicate, "wt", compresslevel=4)
    else:
        dup = open(args.pcr_duplicate, "w")
        
    if len(args.reads_with_improper_umi) != 0:
        if re.search("\.gz|\.gzip", args.reads_with_improper_umi):
            imp = gzip.open(args.reads_with_improper_umi, "wt",
                            compresslevel=4)
        else:
            imp = open(args.reads_with_improper_umi, "w")

    stats = umi.SraRunStats()
    if is_gzipped(fn):
        f = gzip.open(fn, "rt")
        sys.stderr.write("Input is gzipped.\n")
    else:
        f = open(fn)
        
    # The two methods (unique and network-based) are implemented here
    if umi_errors_allowed == 0 and not FORCE_GRAPH:
        # Does not allow erros in UMIs
        # Small RNA insert + umi should be unique; otherwise, it is a duplicate
        insert_umi = {}
        c = 0
        while True:
            c += 1
            if c % 4 == 1:
                r_name = f.readline().strip()
                if not r_name:
                    break
                stats["n_read"] += 1
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
                        print(r_name, file=imp)
                        print(r_seq, file=imp)
                        print(r_info, file=imp)
                        print(r_qual, file=imp)
                else:
                    r_name_proc = get_header_with_umi(r_name_proc, r_bc)
                    stats["n_with_proper_umi"] += 1
                    k = r_seq_proc + "_" + r_bc
                    # print k
                    if k in insert_umi:
                        stats["n_duplicate"] += 1
                        print(r_name_proc, file=dup)
                        print(r_seq_proc, file=dup)
                        print(r_info_proc, file=dup)
                        print(r_qual_proc, file=dup)
                    else:
                        stats["n_non_duplicate"] += 1
                        insert_umi[k] = True
                        print(r_name_proc, file=out)
                        print(r_seq_proc, file=out)
                        print(r_info_proc, file=out)
                        print(r_qual_proc, file=out)

        out.close()
        dup.close()

        if len(args.reads_with_improper_umi) != 0:
            imp.close()
        stats.report()
        # sys.stderr.write("\n")
        # sys.stderr.write("Stats: \n")
        # sys.stderr.write("Total input reads:\t" + str(c/4) + "\n")
        # sys.stderr.write("Reads dropped due to improper UMI:\t" +
        #                  str(stats["n_without_proper_umi"]) + "\n")
        # sys.stderr.write("Final proper read:\t" +
        #                  str(stats["n_with_proper_umi"]) + "\n")
        # sys.stderr.write("\tReads that are duplicates:\t" +
        #                  str(stats["n_duplpicate"]) + "\n")
        # sys.stderr.write("\tReads that are non-duplicates:\t" +
        #                  str(stats["n_non_duplicate"]))

        print2("")
        if verbose:
            pass
    else:
        # Use the network-based method
        # A dictionary of dictionaries: insert2umi[insert]["AAAA"] = umi_count
        insert2umis = {}
        insertumi2read = {}
        # rname2read keeps all reads so that we can figure out which reads are
        # PCR duplicates
        rname2read = {}
        c = 0
        while True:
            c += 1
            if c % 4 == 1:
                r_name = f.readline().strip()
                if not r_name:
                    break
                stats["n_read"] += 1
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
                        print(r_name, file=imp)
                        print(r_seq, file=imp)
                        print(r_info, file=imp)
                        print(r_qual, file=imp)
                else:
                    r_name_proc = get_header_with_umi(r_name_proc, r_bc)
                    stats["n_with_proper_umi"] += 1
                    r = umi.SraRead(r_name_proc, r_seq_proc, 
                                    r_info_proc, r_qual_proc, "r1")
                    insertumi2read[r_seq_proc + r_bc] = r
                    rname2read[r_name_proc] = r
                    # Construct insert2umi
                    if r_seq_proc not in insert2umis:
                        insert2umis[r_seq_proc] = {}
                    if r_bc not in insert2umis[r_seq_proc]:
                        insert2umis[r_seq_proc][r_bc] = 1
                    else:
                        insert2umis[r_seq_proc][r_bc] += 1

        # Now that everything is stored in dictionaries, use the graph
        non_duplicate_rname = {}
        for i in insert2umis:
            umis = insert2umis[i]  # umis is a dict: umi to umi count
            # Build the graph if there are 2 or more reads for the same
            # insert seq
            repr_umis = []
            if (len(umis) > 1):
                G = umi_graph.UmiGraph(umis, max_ed=umi_errors_allowed)
                # If debugging mode is on and the UMI graph identifies more duplicates
                # than the unique method, then print something
                if DEBUG and G.number_true_umi() != len(umis):
                    print(umis)
                    print("Number of true UMIs for this locus: {}".
                          format(G.number_true_umi()))
                    print("Representative UMIs for this UMI graph:")
                    print(G.get_repr_umi())
                    print("Subgraphs of weakly connected components")
                    clusters = G.get_umi_clusters()
                    for c in clusters:
                        print(c)
                    print()
                repr_umis = G.get_repr_umi()
            # No need to build the graph if there is just 1 reads for this
            # insert seq
            else:
                repr_umis = list(umis.keys())

            if len(umis) != len(repr_umis):
                print(umis)

            for u in repr_umis:
                r = insertumi2read[i + u]
                stats["n_non_duplicate"] += 1
                out.write(str(r))
                non_duplicate_rname[r.r_name] = True

        # Print duplicates
        for rname in rname2read:
            if rname not in non_duplicate_rname:
                r = rname2read[rname]
                dup.write(str(r))
                stats["n_duplicate"] += 1

        out.close()
        dup.close()
        stats.report()
        
        
if __name__ == "__main__":
    main()
