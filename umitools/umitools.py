#!/usr/bin/env python3

# The main entrance for the toolkit

import argparse
import os
import sys
import subprocess


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

msg = '''Possible subcommands:
use 'extract_small' to deduplicate a UMI small RNA-seq fastq file;
use 'extract' for UMI RNA-seq to put UMIs into read names;
use 'mark' for marking PCR duplicates given a bam/sam file;
'''.format(sys.argv[0])

parser = argparse.ArgumentParser()

subparsers = parser.add_subparsers(help=msg,
                                   dest='subcommand')
p1 = subparsers.add_parser("extract")
p1 = subparsers.add_parser("extract_small")
p1 = subparsers.add_parser("mark")

args = parser.parse_args()
base = os.path.dirname((os.path.realpath(__file__)))

sub_args = sys.argv[2:]
if len(sub_args) == 0:
    sub_args = ["-h"]

# Extract the subcommand and run the corresponding script
if args.subcommand in ("extract",
                       "reformat_umi_rsq", "reformat_umi"
                       "format_umi_rsq", "format_umi",
                       "extract_rsq", "extract_rna_seq"):
    subprocess.call(["reformat_umi_fastq.py"] + sub_args)

elif args.subcommand in ("extract_small",
                         "reformat_umi_sra", "reformat_umi_small_rna",
                         "reformat_umi_small_rna_seq", "extract_sra",
                         "extract_small_rna", "extract_small_rna_seq",
                         "small", "smallrna", "smallrnaseq"):
    subprocess.call(["reformat_umi_sra_fastq.py"] + sub_args)

elif args.subcommand in ("mark", "mark_duplicates", "mark_duplicate",
                         "rsq_mark_duplicates", "rsq_mark_duplicate"):
    subprocess.call(["umi_mark_duplicates.py"] + sub_args)
    




