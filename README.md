# Description
A toolset for handling sequencing data with unique molecular identifiers (UMIs)

# Installation
To install the toolset, run

`pip3 install umitools  # add --user if you want to install it to your own directory`

Skip to the next section, unless you would like to modify it. You can grab the version on GitHub:

`git clone https://github.com/weng-lab/umitools.git` 

, modify it and then install from the modified source: 

`pip3 install -e /path/to/umitools`

# How to process UMI small RNA-seq data
0. (Skip to the next step if small RNA-seq adapters are already removed.) To process a fastq (`raw.fq.gz`) file from your UMI small RNA-seq data, you can first remove the 3' end small RNA-seq adapter. In this example, I use `fastx_clipper` from the [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/) and the adapter sequence is `TGGAATTCTCGGGTGCCAAGG`:

`zcat raw.fq.gz | fastx_clipper -a TGGAATTCTCGGGTGCCAAGG -l 48 -c -Q33 2> raw.clipped.log | gzip -c - > clipped.fq.gz`

`-l 48` specified the minimum length of the reads after the adapter removal, since I want to make sure all reads are at least 18 nt (18 nt + 15 nt in the 5' UMI + 15 nt in the 3' UMI).

1. To identify UMIs, using our test data:
`wget -O clipped.fq.gz "https://github.com/weng-lab/umitools/raw/master/umitools/testdata/umitools.test.sRNA-seq.fq.gz"`

`umi_reformat_sra_fastq -i clipped.fq.gz -o sra.umi.fq -d sra.dup.fq`

Not sure if your libraries have high-quality UMIs at proper positions? Run the following to see which reads have improper UMIs.

`umi_reformat_sra_fastq -i clipped.fq.gz -o sra.umi.fq -d sra.dup.fq --reads-with-improper-umi sra.improper_umi.fq`

# How to process UMI RNA-seq data

1. Say the read1 and read2 files are `r1.fq.gz` and `r2.fq.gz`. In order to identify reads with proper UMIs and parse out their UMIs, you can run:

`umi_reformat_fastq -l r1.fq.gz -r r2.fq.gz -L r1.fmt.fq.gz -R r2.fmt.fq.gz`

And it will give you some stats for your UMI RNA-seq data.

2. Then you can use your favorite RNA-seq aligner (e.g. STAR) to map these reads to the genome and get a BAM/SAM file (e.g., `fmt.bam`). To mark the reads with PCR duplicates, assuming you want to use 8 threads, simply run

`umi_mark_duplicates -f fmt.bam -p 8`

Reads that are identified as PCR duplicates will have the flag `0x400`. If your downstream analysis (e.g., Picard) can take into consideration this flag, then you are good to go! Otherwise, you can remove those PCR duplicates:

`samtools view -b -h -F 0x400 fmt.deumi.sorted.bam > fmt.deumi.F400.sorted.bam`

You can then feed the bam file without PCR duplicates to your downstream analysis.

# How UMI locators are handled
For UMI RNA-seq, the UMI locator in each read is required to exactly match GGG, TCA, or ATC. You can customize the locator sequence by setting `--umi-locator LOCATOR1,LOCATOR2,LOCATOR3,LOCATOR4` when you run `umi_reformat_fastq`.

For UMI small RNA-seq, the default setting requires that the 5\' UMI locator in each read should match `NNNCGANNNTACNNN` or `NNNATCNNNAGTNNN`, AND 3\' UMI locator should match `NNNGTCNNNTAGNNN` where N's are not required to match and there is at most 1 error across all non-N positions. You can customized the locator sequence for small RNA-seq by setting `--umi-pattern-5` and `--umi-pattern-3`. You can further tweak the number of errors allowed by changing `N_MISMATCH_ALLOWED_IN_UMI_LOCATOR` in the script.

# Other utilities

### umi_simulator.py
A simple in silico PCR simulator for UMI reads.

### find_hot_loci.py
This script can find those "hot" loci, i.e. those loci that produce a huge number of reads and then it outputs a histogram. Optionally, you can include -o option so that it also outputs the corresponding bam records.

# Contact us
Yu Fu (Yu.Fu {at} umassmed.edu)

