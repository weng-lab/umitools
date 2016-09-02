# Description
A toolset for handling sequencing data with unique molecular identifiers (UMIs)

# Installation
`pip install pyfaidx  # add --user if you want to install it to your own directory`

# Usage
### reformat_umi_fastq.py
A script to reformat reads in a UMI fastq file so that the name of each record contains the UMI

### umi_mark_duplicates.py
A pair of FASTQ files are first reformatted using reformat_umi_fastq.py and then is aligned to get the bam file. This script can parse the umi barcode in the name of each read to mark duplicates.

### find_hot_loci.py
This script can find those "hot" loci, i.e. those loci that produce a huge number of reads and then it outputs a histogram. Optionally, you can include -o option so that it also outputs the corresponding bam records

### umi_simulator.py
A simple in silico PCR simulator for UMI reads
