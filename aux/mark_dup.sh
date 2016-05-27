#!/usr/bin/env bash

# Given a bam file that contains barcode info in the header, this pipeline mark duplicates
# using umitools and Picard (optional) and output the bam files with duplicates marked
# It will also output extra bam files that contains only duplicates (f400) and exclude duplicates (F400)
# to facilitate downstream analysis

# input file must have .sorted.bam as the suffix

picard="java -XX:ParallelGCThreads=8 -jar /home/fuy2/src/picard-tools-1.138/picard.jar"

while getopts "hi:Rc:g:o:sD" OPTION; do
    case $OPTION in
	h)usage && exit 0 ;;
	i)input_bam=$OPTARG ;;
	c)CPU=$OPTARG ;;
	D)DEBUG=1 ;; ## DEBUG mode: skipping some steps to make testing easier
	g)GENOME=$OPTARG ;;
    esac
done

if [[ -z ${CPU} ]]; then
    CPU=8
fi

if [[ -z ${GENOME} ]]; then
    GENOME=mm10
fi

if [[ ${GENOME} == 'mm10' ]]; then
    GTF=~/data/shared/mm10/gencode.vM4.corrected_chrom_names.gtf 
fi

fn=$(basename ${input_bam})
prefix=${fn%.bam}
prefix=${prefix%.sorted}
echo ${prefix}

umi_mark_duplicates.py -p ${CPU} -f ${input_bam} 2>${i}.deumi.log
deumi_bam=${prefix}.deumi.sorted.bam

# F400 and f400 bam files based on UMI info
if [[ ! -f ${prefix}.deumi.ok ]]; then
    samtools view -b -h -F 0x400 ${deumi_bam} > ${prefix}.F400.sorted.bam && \
	samtools view -b -h -f 0x400 ${deumi_bam} > ${prefix}.f400.sorted.bam && \
	touch ${prefix}.deumi.ok
else
    echo "Skipping the deumi step. If you want to redo it, please delete ${prefix}.deumi.ok."
fi

# Picard
if [[ ! -f ${prefix}.picard.ok ]]; then
    # F400 and f400 bam files based on Picard
    ${picard} MarkDuplicates I=${prefix}.sorted.bam O=${prefix}.picard_dup_marked.bam M=${prefix}.picard_dup_marked.metric &> ${prefix}.picard_dup_marked.log && \
	samtools view -b -F 0x400 ${prefix}.picard_dup_marked.bam > ${prefix}.picard.F400.bam && \
	touch ${prefix}.picard.ok
else
    echo "Skipping Picard step. If you want otredo it, please delete ${prefix}.picard.ok."
fi

# Counting
cmds=jobs.${RANDOM}${RANDOM}.sh

if [[ ! -f ${prefix}.count.ok ]]; then
    echo "htseq_bam_to_count.sh ${GTF} ${prefix}.F400.sorted.bam &>${i}.log" > ${cmds}
    echo "htseq_bam_to_count.sh ${GTF} ${prefix}.f400.sorted.bam &>${i}.log" > ${cmds}
    echo "htseq_bam_to_count.sh ${GTF} ${prefix}.picard.F400.bam &>${i}.log" > ${cmds}
    cat ${cmds} | parallel -j ${CPU} --progress --dry-run && touch ${prefix}.count.ok
else
    echo "Skipping counting step. If you want to redo it, please delete ${prefix}.count.ok"
fi

echo "Done"
