#!/usr/bin/env python3

# The main entrance for the toolkit

import argparse
import os
import sys
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('subcommand', help="foo help")

args = parser.parse_args()
base = os.path.dirname((os.path.realpath(__file__)))

# Extract the subcommand and run the corresponding script
if args.subcommand in ("reformat_umi_rsq", "reformat_umi"
                       "format_umi_rsq", "format_umi"):
    sub_args = sys.argv[2:]
    subprocess.call(["reformat_umi_fastq.py", ] + sub_args)

elif args.subcommand in ("reformat_umi_sra", "reformat_umi_small_rna",
                         "reformat_umi_small_rna_seq"):
    sub_args = sys.argv[2:]

else:
    sys.stderr.write("Unknown subcommand")
    parser.print_help()
