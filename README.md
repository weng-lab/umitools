## Description
A toolset for handling sequencing data with unique molecular identifiers (UMIs)

## Installation
This tools set requires Python 3.

To install `umitools`, run

```shell
pip3 install umitools  # add --user if you want to install it to your own directory
```

## How to process UMI small RNA-seq data
#### 0. (Skip to the next step if you have data.) Download the test data

```shell
wget -O clipped.fq.gz "https://github.com/weng-lab/umitools/raw/master/umitools/testdata/umitools.test.sRNA-seq.fq.gz"
```

#### 1. Identify UMIs:

```shell
umi_reformat_sra_fastq -i clipped.fq.gz -o sra.umi.fq -d sra.dup.fq
```


## How to process UMI RNA-seq data
#### 0. (Skip to the next step if you have data.) Download the test data

```shell
wget -O "r1.fq.gz" "https://github.com/weng-lab/umitools/raw/master/umitools/testdata/umitools.test.RNA-seq.r1.fq.gz"
wget -O "r2.fq.gz" "https://github.com/weng-lab/umitools/raw/master/umitools/testdata/umitools.test.RNA-seq.r2.fq.gz"
```

#### 1. To identify reads with proper UMIs and parse out their UMIs, you can run:

```shell
umi_reformat_fastq -l r1.fq.gz -r r2.fq.gz -L r1.fmt.fq.gz -R r2.fmt.fq.gz
```

And it will output some stats for your UMI RNA-seq data.

#### 2. Then you can use your favorite RNA-seq aligner (e.g. STAR) to map these reads to the genome and get a BAM/SAM file (e.g., `fmt.bam`). To download an example, run

```shell
wget -O fmt.bam https://github.com/weng-lab/umitools/raw/master/umitools/testdata/umitools.test.RNA-seq.sorted.bam
```

To mark the reads with PCR duplicates (and assuming you want to use 8 threads), run

```shell
umi_mark_duplicates -f fmt.bam -p 8
```

And it will produce `fmt.deumi.sorted.bam` in which reads that are identified as PCR duplicates will have the flag `0x400`. If your downstream analysis (e.g., Picard) can take into consideration this flag, then you are good to go! Otherwise, you can just eliminate PCR duplicates:

```shell
samtools view -b -h -F 0x400 fmt.deumi.sorted.bam > fmt.deumi.F400.sorted.bam
```

You can then feed the bam file without PCR duplicates to your downstream analysis.

## How UMI locators are handled
For UMI RNA-seq, the UMI locator in each read is required to exactly match GGG, TCA, or ATC. You can customize the locator sequence by setting `--umi-locator LOCATOR1,LOCATOR2,LOCATOR3,LOCATOR4` when you run `umi_reformat_fastq`.

For UMI small RNA-seq, the default setting requires that the 5\' UMI locator in each read should match `NNNCGANNNTACNNN` or `NNNATCNNNAGTNNN`, AND 3\' UMI locator should match `NNNGTCNNNTAGNNN` where N's are not required to match and there is at most 1 error across all non-N positions. You can customized the locator sequence for small RNA-seq by setting `--umi-pattern-5` and `--umi-pattern-3`. You can further tweak the number of errors allowed by changing `N_MISMATCH_ALLOWED_IN_UMI_LOCATOR` in the script.

## Other utilities

#### umi_simulator.py
A simple in silico PCR simulator for UMI reads.

#### find_hot_loci.py
This script can find those "hot" loci, i.e. those loci that produce a huge number of reads and then it outputs a histogram. Optionally, you can include -o option so that it also outputs the corresponding bam records.

## FAQ 

#### How to remove 3' end small RNA-seq adapter
There are many tools to remove adapters. This is just one example. To process a fastq (`raw.fq.gz`) file from your UMI small RNA-seq data, you can first remove the 3' end small RNA-seq adapter. For example, you can use `fastx_clipper` from the [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/) and the adapter sequence is `TGGAATTCTCGGGTGCCAAGG`:

```shell
zcat raw.fq.gz | fastx_clipper -a TGGAATTCTCGGGTGCCAAGG -l 48 -c -Q33 2> raw.clipped.log | gzip -c - > clipped.fq.gz
```

where `-l 48` specified the minimum length of the reads after the adapter removal, since I want to make sure all reads are at least 18 nt (18 nt + 15 nt in the 5' UMI + 15 nt in the 3' UMI).


### Not sure if your libraries have high-quality UMIs at proper positions? 

To see which reads have improper UMIs, run

```shell
umi_reformat_sra_fastq -i clipped.fq.gz -o sra.umi.fq -d sra.dup.fq --reads-with-improper-umi sra.improper_umi.fq
```
where `sra.umi.fq` contains all the non-duplicate reads and `sra.dup.fq` contains all duplicates.

#### Feeling adventurous? You can install the git version
1. Grab the version on GitHub:

```shell
git clone https://github.com/weng-lab/umitools.git
```

2. Install it in editable mode: 

```shell
pip3 install -e /path/to/umitools
```

## Citation
Fu, Y., Wu, P., Beane, T., Zamore, P.D., and Weng, Z. (2018). Elimination of PCR duplicates in RNA-seq and small RNA-seq using unique molecular identifiers. BioRxiv 251892.

## Contact us
Yu Fu (Yu.Fu {at} umassmed.edu)

