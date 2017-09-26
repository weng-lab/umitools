#!/usr/bin/env python3

# Estimate the error rates of short reads using the fixed portion of reads
# in UMI small RNA-seq

# TODO: Do I need to do this for UMI RNA-seq?

import umi
import argparse



desc = "A script that estimates the sequencing error rates"
parser = argparse.ArgumentParser(description=desc,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

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

ui = umi.SraUmiInfo(umi_pat5, umi_pat3)



