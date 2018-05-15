#!/usr/bin/env python3

# TODO: The main entrance for the toolkit

import argparse
import os
import sys
import subprocess


def main():
    parser = argparse.ArgumentParser()
    msg = '''
    For UMI RNA-seq:
    First, use umitools reformat_fastq to identify UMIs in UMI RNA-seq
    Then, use umitools umi_mark_duplicates to mark PCR duplicates

    For UMI small RNA-seq:
    Use umitools reformat_sra_fastq to identify UMIs and PCR duplicates

    To simulate UMIs, use umitools simulate.
    '''.format(sys.argv[0])

    parser.add_argument(help=msg,
                        dest='subcommand')

    # base = os.path.dirname((os.path.realpath(__file__)))

    # If there is not parameters, just print the help and exit
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # If there are parameters, consider first one as the subcommand
    d = vars(parser.parse_args(sys.argv[1:2]))
    if 'subcommand' in d:
        subcmd = d['subcommand']

    print(subcmd)

    if len(sys.argv) > 2:
        params = sys.argv[2:]

    else:
        params = []

    # Even though 'reformat_fastq', 'mark_duplicates' and 'reformat_sra_fastq'
    # are the advertised subcomamnds, it is also able to handle similar subcommands
    if subcmd in ["reformat_fastq", "umi_reformat_fastq"]:
        subprocess.call(["umi_reformat_fastq"] + params)
    elif subcmd in ["mark_duplicates"]:
        subprocess.call(["umi_mark_duplicates"] + params)
    elif subcmd in ["reformat_sra_fastq", "umi_reformat_sra_fastq"]:
        subprocess.call(["umi_reformat_sra_fastq"] + params)
    else:
        sys.stderr.write("Unknown commands!")
        parser.print_help()
        sys.exit(1)
