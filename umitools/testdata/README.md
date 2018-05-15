
# umitools.test.RNA-seq.log contains the correct log from 

`umi_reformat_fastq -l umitools.test.RNA-seq.r1.fq.gz -r umitools.test.RNA-seq.r2.fq.gz -L /tmp/1.fastq.gz -R /tmp/2.fastq.gz`

# umitools.test.RNA-seq.sorted.bam.stats contains stats after PCR duplicate identification

`umi_mark_duplicates -f umitools.test.RNA-seq.sorted.bam`
`samtools flagstat umitools.test.RNA-seq.sorted.deumi.sorted.bam`


# umitools.test.sRNA-seq.log contains the correct log from

`umi_reformat_sra_fastq -i umitools.test.sRNA-seq.fq.gz -o /tmp/output.fq.gz -d /tmp/dup.fq.gz`
