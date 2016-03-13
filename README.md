# umitools
A toolset for handling sequencing data with unique molecular identifiers (UMIs)

### reformat_umi_fastq.py
A script to reformat reads in a UMI fastq file so that the name of each record contains the UMI

### umi_mark_duplicates.py
A pair of FASTQ files are first reformatted using reformat_umi_fastq.py and then is aligned to get the bam file. This script can parse the umi barcode in the name of each read to mark duplicates.

find_hot_loci.py

umi_loci_with_duplicates.py

umi_simulator.py
