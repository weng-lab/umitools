#!/usr/bin/env python3

# The main entrance for the toolkit

import argparse
import os
import sys
import subprocess


parser = argparse.ArgumentParser()
msg = '''Use 'extract_small' to deduplicate a UMI small RNA-seq fastq file;
use 'extract' for UMI RNA-seq to put UMIs into read names;
use 'mark' for marking PCR duplicates given a bam/sam file;
'''.format(sys.argv[0])

subparsers = parser.add_subparsers(help=msg,
                                   dest='subcommand')
p1 = subparsers.add_parser("extract")
# a1 = ("extract", "reformat_umi_rsq", "reformat_umi"
#       "format_umi_rsq", "format_umi",
#       "extract_rsq", "extract_rna_seq")

p2 = subparsers.add_parser("extract_small")
# a2 = ("extract_small",
#       "reformat_umi_sra", "reformat_umi_small_rna",
#       "reformat_umi_small_rna_seq", "extract_sra",
#       "extract_small_rna", "extract_small_rna_seq",
#       "small", "smallrna", "smallrnaseq")

p3 = subparsers.add_parser("mark")
# a3 = ("mark", "mark_duplicates", "mark_duplicate",
#       "rsq_mark_duplicates", "rsq_mark_duplicate")

args = parser.parse_args()
base = os.path.dirname((os.path.realpath(__file__)))

if len(sys.argv) == 1:
    parser.print_help()

sub_args = sys.argv[2:]
if len(sub_args) == 0:
    sub_args = ["-h"]

# This allows the variants of options to work
if args.subcommand == "extract":
    subprocess.call(["reformat_umi_fastq.py"] + sub_args)

elif args.subcommand == "extract_small":
    subprocess.call(["reformat_umi_sra_fastq.py"] + sub_args)

elif args.subcommand == "mark":
    subprocess.call(["umi_mark_duplicates.py"] + sub_args)
