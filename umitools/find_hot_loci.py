#!/usr/bin/env python

import argparse
import sys
import pysam

global DEBUG
DEBUG = False


def read_to_key(read):
    # read_chr_id = read.reference_id
    read_chr = read.reference_name
    read5 = read.reference_start
    read_strand = ''
    if read.is_reverse:
        read_strand = '-'
    else:
        read_strand = '+'
    # We need abs() here as we do not know if r1 is the leftmost read or not
    tl = abs(read.template_length)
    k = ''.join((str(read_chr), ":", str(read5), "-", str(tl), read_strand))
    # print k
    return k


def get_reads_for_each_pos(fn):
    c = 0
    bam = pysam.AlignmentFile(fn, "rb")
    pos2read = {}
    for read in bam.fetch():
        c += 1
        if c % 100000 == 0:
            print >>sys.stderr, "Processed " + str(c) + " entries..."
        if read.is_read1:
            read_n = read.query_name
            # read_bc = read_n.split("_")[1]
            read_name = read_n.split("_")[0]
            k = read_to_key(read)
            if k in pos2read:
                pos2read[k].append(read_name)
            else:
                pos2read[k] = [read_name, ]
    return pos2read


def output_bam(fn, pos2read, cutoff=100):
    '''For those loci that generates lots of reads, this function returns
    the corresponding records'''
    all_read_names = {}
    for i in pos2read:
        if len(pos2read[i]) >= cutoff:
            for j in pos2read[i]:
                all_read_names[j] = 0
    bam = pysam.AlignmentFile(fn, "rb")
    out = pysam.AlignmentFile(fn + ".hot_loci.bam", "wb", template=bam)
    c = 0
    for read in bam.fetch():
        c += 1
        if c % 100000 == 0:
            print >>sys.stderr, "Processed " + str(c) + " entries..."
        read_n = read.query_name
        # read_bc = read_n.split("_")[1]
        read_name = read_n.split("_")[0]
        if read_name in all_read_names:
            out.write(read)
    out.close()
            
    
def print_hist(pos2read):
    hist = {}
    for k in pos2read:
        l = len(pos2read[k])
        if l in hist:
            hist[l] += 1
        else:
            hist[l] = 1
    for i in range(max(hist.keys()) + 1):
        if i in hist:
            print str(i) + "\t" + str(hist[i])
        else:
            print str(i) + "\t" + str(0)

            
def get_hottest_n(pos2read, n):
    '''Returns a list of hot loci (number of reads and position)
'''
    ret = []
    lens = [len(pos2read[i]) for i in pos2read]
    tmp = sorted(lens, reverse=True)
    my_min = min(tmp[1:100])
    for i in pos2read:
        if len(pos2read[i]) >= my_min:
            ret.append(str(len(pos2read[i])) + "\t" + i)
    return(ret)


def main():
    parser = argparse.ArgumentParser(description='This script can \
    find those "hot" loci, i.e. those loci that produce a huge number \
    of reads and then it outputs a histogram. Optionally, you can include \
    -o option so that it also outputs the corresponding bam records')
    parser.add_argument('-f', '--file', help='the input bam file',
                        required=True)
    parser.add_argument('-d', '--debug', help='turn on debug mode',
                        action="store_true")
    parser.add_argument('-o', '--output-bam', help='output the bam file \
    containing reads of interest. The default is to output the bam records \
    with more than 100 reads',
                        action="store_true")
    parser.add_argument('--output-bam-cutoff', help='output the bam file \
    containing reads of interest. The default is to output the bam records \
    with more than 100 reads', type=int,
                        default=100)

    parser.add_argument('--output-hottest', help='output the N hottest \
    positions',
                        type=str,
                        default=None)

    parser.add_argument('--hottest-n', help='specify how many hot loci \
    are output. By default, it outputs the hottest 100 loci',
                        type=int,
                        default=100)

    args = parser.parse_args()
    fn = args.file
    
    output_hottest = args.output_hottest
    n_hottest = args.hottest_n
    
    # fn = "Zamore.RSQ.UMI.B6.brain.17dpp.fmt.x_rRNA.mm10g.F400.sorted.bam"
    output_bam_flag = False
    output_bam_flag = args.output_bam
    output_bam_cutoff = args.output_bam_cutoff
    bam_tmp = pysam.AlignmentFile(fn, "rb")
    if not bam_tmp.has_index():
        print >>sys.stderr, "Input file %s is being indexed." % fn
        pysam.index(fn)
    global DEBUG
    DEBUG = args.debug
    pos2read = get_reads_for_each_pos(fn)
    print >>sys.stderr, "The cutoff for hot loci: %d" % output_bam_cutoff
    if output_bam_flag:
        # output_bam(fn, pos2read, cutoff=output_bam_cutoff)
        output_bam(fn, pos2read, cutoff=output_bam_cutoff)
    print_hist(pos2read)

    if output_hottest:
        fh_hot = open(output_hottest, "w")
        print >>fh_hot, "\n".join(get_hottest_n(pos2read, n_hottest))

if __name__ == '__main__':
    main()
