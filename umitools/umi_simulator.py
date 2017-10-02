#!/usr/bin/env python3

import random
import sys
import argparse

__author__ = "Yu Fu"
__license__ = "GPLv3"

# This is where we define the exp level of each locus
# for the sake of completeness, I also include a locus that has 0 expression
# Note that no two loci can share the same exp level
# LOCI_EXP = range(0, 1000) * 30 + [1000, 3000, 5000, 10000, 30000, 50000, 100000]

# Final one to use
# LOCI_EXP = range(0, 2000) + range(2000, 100000, 500)
LOCI_EXP = list(range(0, 2000)) + list(range(2000, 100000, 500))

# Testing
LOCI_EXP = list(range(0, 100))

# Limit the pool size: if after any PCR cycle, the pool is larger that this, then downsample the pool
# This should be larger than final_pool_size
# MAX_POOL_SIZE = 10000000
MAX_POOL_SIZE = 10000000


def print2(a):
   print(a, file=sys.stderr)

    
class UMIRead:
    '''A class representing UMI reads.
Species represents the ground truth and seq represents the actual sequence
which may or may not have mutations
gid is the true number of molecules initially. species is the true UMI
initially. seq is the sequenced UMI, which may have errors caused by
PCR or sequencing. gid is the ID and also the number of molecules before
PCR amplification.
'''
    def __init__(self, *args):
        self.gid = args[0]
        self.species = args[1]
        if len(args) == 3:
            self.seq = args[2]
        else:
            self.seq = args[1]
        
    def amplify(self, prob=0.8, error=1e-4):
        ## Default PCR amplification efficiency
        mut = {'A': ('C', 'G', 'T'),
               'C': ('A', 'G', 'T'),
               'G': ('A', 'C', 'T'),
               'T': ('A', 'C', 'G')
               }
        # The makes the PCR amplification probability random
        cur_prob = random.uniform(prob, 1)
        if random.random() <= cur_prob:
            if error == 0:
                return (UMIRead(self.gid, self.species, self.seq), 
                        UMIRead(self.gid, self.species, self.seq))
            else:
                new = []
                for i in range(len(self.seq)):
                    if random.random() < error:
                        tmp = random.choice(mut[self.seq[i]])
                        new.append(tmp)
                    else:
                        new.append(self.seq[i])
                # The species after amplification never changes
                return (UMIRead(self.gid, self.species, self.seq), 
                        UMIRead(self.gid, self.species, ''.join(new)))
        else:
            return (UMIRead(self.gid, self.species, self.seq), )

    def __str__(self):
        '''Print the real sequence (before PCR) and the 
final sequence (after PCR and sequencing)
'''
        return "{}\t{}".format(self.species, self.seq)

    def __repr__(self):
        return "UMI: id %s species %s seq %s" % \
            (self.gid, self.species, self.seq)

    
def add_sequencing_error(r, error=0.01):
    '''Returns a new UMIRead object containing the sequencing error
'''
    mut = {'A': ('C', 'G', 'T'),
    'C': ('A', 'G', 'T'),
    'G': ('A', 'C', 'T'),
    'T': ('A', 'C', 'G')
           }
    if error == 0:
        return (UMIRead(r.gid, r.species, r.seq))
    else:
        new = []
        for i in range(len(r.seq)):
            if random.random() < error:
                tmp = random.choice(mut[r.seq[i]])
                new.append(tmp)
            else:
                new.append(r.seq[i])
        # Returns a new UMIRead object. The species after amplification never changes
        return (UMIRead(r.gid, r.species, ''.join(new)))


def test4():
    print("Testing sequencing error (1%, which makes things obvious)")
    c = 0
    n = 10000
    for i in range(n):
        umi = UMIRead("0", 'ACGT')
        new = add_sequencing_error(umi, 0.01)
        if new.seq != umi.seq:
            c += 1
    print("Number of reads out ouf %d with errors: %d" % (n, c))
    print("It should be around 400")

    
def test1():
    ## PCR eff is 1, error rate is 1%
    print("PCR eff is 1, error rate is 1%")
    umi = UMIRead("0", 'ACGT')
    myerr = [0,] * len(umi.species)
    for i in range(10000):
        tmp = umi.amplify(prob=1, error=0.01)
        if len(tmp) > 1:
            my = tmp[1]
            ## print my, umi.species
            for i in range(len(my.seq)):
                if my.seq[i] != umi.species[i]:
                    myerr[i] += 1
    for i in range(len(umi.species)):
        print(str(i) + "\t" + str(myerr[i]))


def test2():
    ## PCR eff is 0.5, error rate is 1%
    ## test if the number of reads doubles
    print("PCR eff is 0.5, error rate is 1%")
    n = 10000
    umi = UMIRead("0", 'ACGT')
    myerr = [0,] * len(umi.species)
    for i in range(10000):
        tmp = umi.amplify(prob=0.5, error=0.01)
        if len(tmp) > 1:
            n += 1
            my = tmp[1]
            ## print my, umi.species
            for i in range(len(my.seq)):
                if my.seq[i] != umi.species[i]:
                    myerr[i] += 1
    for i in range(len(umi.species)):
        print(str(i) + "\t" + str(myerr[i]))
    print("Number of reads in the end: %d" % n)

    
def test3():
    print(UMIRead("0", "ATCG", "ATTT"))


def main():
    parser = argparse.ArgumentParser(description='A simple in silico PCR simulator. It creates an initial set of molecules with absolute numbers ranging from 1 to 100000), simulates PCR and sequencing and and outputs the stats.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', '--pcr-cycle', help='number of PCR cycles', required=True, type=int)
    parser.add_argument('-l', '--umi-length', help='length of UMI', required=True, type=int)
    parser.add_argument('--reads-single-locus', help='number of reads for simulating one locus. Using this option causes the scripts to simulate just one locus', required=False, type=int)
    ## parser.add_argument('-s', '--pool-size', help='initial pool size (number of molecules before PCR)', required=True, type=int)
    ## parser.add_argument('-o', '--output-size', help='final pool size (sequencing depth)', required=True, type=int)
    parser.add_argument('-a', '--amplification-rate', help='successful rate of PCR amplification', required=False, type=float, default=0.7)
    parser.add_argument('--pcr-error', help='error rate of PCR amplification', required=False, type=float, default=3e-5)
    parser.add_argument('--sequencing-error', help='error rate of sequencing', required=False, type=float, default=0.001)    
    args = parser.parse_args()
    
    pool = []
    # k = 4
    k = args.umi_length
    # pool_size = 100
    # pool_size = args.pool_size
    global LOCI_EXP
    if args.reads_single_locus:
        LOCI_EXP = [args.reads_single_locus]
        print2("Only one locus will be simulated!")
    # final_pool_size = 1000
    # final_pool_size = args.output_size
    final_pool_size = 10000000
    # pcr_n = 10
    pcr_n = args.pcr_cycle
    success_rate = args.amplification_rate
    pcr_error = args.pcr_error
    sequencing_error = args.sequencing_error
    # for i in range(pool_size):
    #     tmp = ''.join(random.choice(('A', 'C', 'G', 'T')) for _ in range(k))
    #     pool.append(UMIRead("0", tmp))
    print2("Total number of molecules to start with %d" % sum(LOCI_EXP))
    for i in range(len(LOCI_EXP)):
        for j in range(LOCI_EXP[i]):
            tmp = ''.join(random.choice(('A', 'C', 'G', 'T')) for _ in range(k))
            pool.append(UMIRead(LOCI_EXP[i], tmp))
    for i in range(pcr_n):
        new_pool = []
        for p in pool:
            new_pool.extend(p.amplify(success_rate, pcr_error))
        print2("Done PCR cycle %d." % (i+1))
        print2("Current number of molecules %d." % (len(new_pool)))
        if len(new_pool) <= MAX_POOL_SIZE:
            pool = new_pool
        else:
            pool = random.sample(new_pool, MAX_POOL_SIZE)
    print2("Pool size after PCR: %d" % len(pool))
    if len(pool) > final_pool_size:
        pool2 = random.sample(pool, final_pool_size)
    else:
        pool2 = pool
        print2('''If you are NOT debugging this script, \
chances are some parameters are NOT set correctly.''')
    print2("Number of reads to be sequenced: %d" % len(pool2))
    n_err_read = 0
    for i in pool2:
        if i.species != i.seq:
            n_err_read += 1
    print2("Number of reads with error(s) after PCR: %d" % n_err_read)
    final_pool = []
    for i in pool2:
        r = add_sequencing_error(i, sequencing_error)
        final_pool.append(r)
    print2("Number of reads with error(s) after sequencing: %d" % n_err_reads(final_pool))

    ## Print the stats of the original reads for each loci and final reads for each loci
    loci_reads = {}
    # print LOCI_EXP
    for i in range(len(LOCI_EXP)):
        loci_reads[LOCI_EXP[i]] = []
    for p in final_pool:
        loci_reads[p.gid].append(p)
    # print("%s\t%s\t%s" % ("initial_true_n_fragments", "final_n_reads", 
    #                       "final_n_reads_unique", "final_n_reads_1err", 
    #                       "final_n_reads_2err"))

    ########################################################################
    # Allow errors in UMIs. Graph method
    # TODO: test this block
    ########################################################################
    locus2umi_count_net = {}
    locus2umis = {}
    for i in LOCI_EXP:
        locus2umis[i] = []
        locus2umi_count_net[i] = 0
    for i in final_pool:
        locus = i.gid
        umi = i.seq
        locus2umis[locus].append(umi)
    for locus in locus2umis:
        if len(umis) == 1:
            locus2umi_count_net[locus] = 1
        else:
            umis = locus2umis[locus]
            G = umi_graph.UmiGraph(umis, max_ed=1)
            repr_umis = G.get_repr_umi()
            locus2umi_count_net[locus] = len(repr_umis)
   ########################################################################

        
    for i in range(len(LOCI_EXP)):
        j = loci_reads[LOCI_EXP[i]]
        n_umi = len(set([jj.species for jj in j]))
        n_reads = len(j)
        print("%d\t%d\t%d" % (LOCI_EXP[i], n_reads, n_umi))


    for gid in LOCI_EXP:
        print(i)

    # DO NOT use this
    # Print all reads (The output is impossiblely huge when I need to do thousands of
    # simulations
    # print("#true_seq\tread")
    # for i in final_pool:
    #     print(i)
        

def n_err_reads(p):
    n = 0
    for i in p:
        if i.species != i.seq:
            n += 1
    return n
    

if __name__ == "__main__":
    # test2()
    main()
