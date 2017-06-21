# Description
A toolset for handling sequencing data with unique molecular identifiers (UMIs)

# Installation
`pip install pyfaidx  # add --user if you want to install it to your own directory`

If you would like to modify it, simply grab the version on GitHub:

`git clone https://github.com/weng-lab/umitools.git` and use the python scripts in umitools/umitools.

# How to process UMI small RNA-seq data
To process a fastq (`raw.fq.gz`) file from your UMI small RNA-seq data, you can first remove the 3' end small RNA-seq adapter. In this example, I use fastx_clipper from the FASTX-Toolkit and the adapter sequence is `TGGAATTCTCGGGTGCCAAGG`:

`zcat raw.fq.gz | fastx_clipper -a TGGAATTCTCGGGTGCCAAGG -l 21 -c -Q33 2> raw.clipped.log | gzip -c - > clipped.fq.gz`

To identify UMIs, you can run

`reformat_umi_sra_fastq.py -i clipped.fq.gz -o sra.umi.fq -d sra.dup.fq`

Not sure if your libraries have high-quality UMIs at proper positions? Run the following to see which reads have improper UMIs.

`reformat_umi_sra_fastq.py -i clipped.fq.gz -o sra.umi.fq -d sra.dup.fq --reads-with-improper-umi sra.improper_umi.fq`

# Scripts for UMI RNA-seq
### reformat_umi_fastq.py
A script to reformat reads in a UMI fastq file so that the name of each record contains the UMI.

### umi_mark_duplicates.py
A pair of FASTQ files are first reformatted using reformat_umi_fastq.py and then is aligned to get the bam file. This script can parse the umi barcode in the name of each read to mark duplicates.

### find_hot_loci.py
This script can find those "hot" loci, i.e. those loci that produce a huge number of reads and then it outputs a histogram. Optionally, you can include -o option so that it also outputs the corresponding bam records.

### umi_simulator.py
A simple in silico PCR simulator for UMI reads.

# Scripts for UMI small RNA-seq
### reformat_umi_sra_fastq.py
This script idenitifies UMIs from UMI small RNA-seq data.
