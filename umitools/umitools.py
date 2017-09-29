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

parser.add_argument(help=msg,
                    dest='subcommand')

subcommand, params = parser.parse_known_args()

a1 = ("extract", "reformat_umi_rsq", "reformat_umi"
      "format_umi_rsq", "format_umi",
      "extract_rsq", "extract_rna_seq")

a2 = ("extract_small",
      "reformat_umi_sra", "reformat_umi_small_rna",
      "reformat_umi_small_rna_seq", "extract_sra",
      "extract_small_rna", "extract_small_rna_seq",
      "small", "smallrna", "smallrnaseq")

a3 = ("mark", "mark_duplicates", "mark_duplicate",
      "rsq_mark_duplicates", "rsq_mark_duplicate")

args = sys.argv
base = os.path.dirname((os.path.realpath(__file__)))

if len(sys.argv) == 1:
    print(msg)
    sys.exit(0)

sub_args = sys.argv[2:]
if len(sub_args) == 0:
    sub_args = ["-h"]

subcommand = args[1]

# This allows the variants of options to work
if subcommand in a1:
    subprocess.call(["reformat_umi_fastq.py"] + sub_args)

elif subcommand in a2:
    subprocess.call(["reformat_umi_sra_fastq.py"] + sub_args)

elif subcommand in a3:
    subprocess.call(["umi_mark_duplicates.py"] + sub_args)
