#!/usr/bin/env python3

import random
import sys
import argparse
import umitools.umi_graph as umi_graph
import multiprocessing as mp

__author__ = "Yu Fu"
__license__ = "GPLv3"

# This is where we define the exp level of each locus
# for the sake of completeness, I also include a locus that has 0 expression
# Note that no two loci can share the same exp level
# LOCI_EXP = range(0, 1000) * 30 + [1000, 3000, 5000, 10000, 30000, 50000, 100000]

# Final one to use
# LOCI_EXP = range(0, 2000) + range(2000, 100000, 500)
# LOCI_EXP = list(range(1, 2000)) + list(range(2000, 100000, 500))

# Testing
# LOCI_EXP = list(range(1, 100))

# Limit the pool size: if after any PCR cycle, the pool is larger that this, then downsample the pool
# This should be larger than final_pool_size
global MAX_POOL_SIZE
MAX_POOL_SIZE = 1000000


def n_err_reads(p):
    n = 0
    for i in p:
        if i.species != i.seq:
            n += 1
    return n
    

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
        
    def amplify(self, prob=0.8, error=3e-5):
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

    
def add_sequencing_error(r, error=0.001):
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
    parser = argparse.ArgumentParser(description='A simple in silico PCR simulator. It creates an initial set of molecules for one locus, simulates PCR and sequencing and outputs the stats. It is necessary to specify the seed, since by default, this simulator uses 0 as the seed', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', '--pcr-cycle', help='number of PCR cycles', required=False, default=10, type=int)
    parser.add_argument('-l', '--umi-length', help='length of UMI', required=False, default=18, type=int)
    parser.add_argument('-s', '--pool-size', help='initial pool size (number of molecules before PCR)', required=False, type=int, default=100)
    parser.add_argument('-o', '--output-size', help='final pool size (sequencing depth, i.e. number of reads sampled from the PCR amplified pool)', required=False, type=int, default=100)
    parser.add_argument('-a', '--amplification-rate', help='successful rate of PCR amplification. The actual amplification rate is uniformally distributed between this number and 1', required=False, type=float, default=0.8)
    parser.add_argument('--pcr-error', help='error rate of PCR amplification', required=False, type=float, default=3e-5)
    parser.add_argument('--sequencing-error', help='error rate of sequencing', required=False, type=float, default=0.001)
    parser.add_argument('--task', help='tasks to simulate multiple conditions. Other PCR-related arguments are ignored when this option is used', required=False, type=str)
    parser.add_argument('--task-rep', help='number of replicates when doing tasks', required=False, type=int, default=100)    
    parser.add_argument('--cpu', help='Tasks supporting multiprocessing.', required=False, default=16)        
    parser.add_argument('--seed', help='starting seed', type=int, required=False, default=0)
    # parser.add_argument('--reads-single-locus', help='number of reads for simulating one locus. Using this option causes the scripts to simulate just one locus', required=False, type=int)
            
    args = parser.parse_args()
    # k = 4
    k = args.umi_length
    pool_size = args.pool_size
    final_pool_size = args.output_size
    pcr_n = args.pcr_cycle
    success_rate = args.amplification_rate
    pcr_error = args.pcr_error
    sequencing_error = args.sequencing_error
    n_cpu = int(args.cpu)
    task = args.task
    task_rep = args.task_rep
    seeds = range(args.seed, args.seed + task_rep)
    
    print("# Input conditions:")
    print("# UMI length: {}".format(k))
    print("# PCR cycle: {}".format(pcr_n))
    print("# Initial pool size: {}".format(pool_size))
    print("# Final pool size: {}".format(final_pool_size))
    print("# PCR success rate: between {} and 1".format(success_rate))
    print("# Sequencing error rate: {}".format(sequencing_error))
    if args.task:
        print("# Task: {}".format(args.task))
    else:
        print("# Task: {}".format("just one run"))
    print("# Number for replicates for this task: {}".format(args.task_rep))
    print("# Seeds: ", end="")
    print(list(seeds))
    
    if task is None:
        # simulate(pool_size, final_pool_size, k, pcr_n,
        #          success_rate, pcr_error, sequencing_error, args.seed)
        print2("No task specified. All parameters are set as specified.")
        simulate_multiple(pool_size, final_pool_size, k, pcr_n,
                          success_rate, pcr_error, sequencing_error,
                          n_cpu=n_cpu, seeds=seeds)
        
    elif task == "pcr_cycle":
        print2("Variable PCR cycles. Other parameters are set as specified.")
        for pcr_n in range(1, 21):
            simulate_multiple(pool_size, final_pool_size, k, pcr_n,
                              success_rate, pcr_error, sequencing_error,
                              n_cpu=n_cpu, seeds=seeds)

    elif task == "pcr_cycle_high":
        print2("Variable PCR cycles. Other parameters are set as specified.")
        for pcr_n in range(21, 31):
            simulate_multiple(pool_size, final_pool_size, k, pcr_n,
                              success_rate, pcr_error, sequencing_error,
                              n_cpu=n_cpu, seeds=seeds)
            
    elif task == "umi_length":
        for k in range(4, 22):
            simulate_multiple(pool_size, final_pool_size, k, pcr_n,
                              success_rate, pcr_error, sequencing_error,
                              n_cpu=n_cpu, seeds=seeds)

    elif task == "pcr_err":
        a = [x / 10.0 for x in range(-70, -29, 1)]
        a = [10**i for i in a]
        for pcr_error in a:
            simulate_multiple(pool_size, final_pool_size, k, pcr_n,
                              success_rate, pcr_error, sequencing_error,
                              n_cpu=n_cpu, seeds=seeds)
    
    elif task == "sequencing_error":
        a = [i / 10 for i in range(-50, -9, 1)]
        a = [10 ** i for i in a]
        for sequencing_error in a:
            simulate_multiple(pool_size, final_pool_size, k, pcr_n,
                              success_rate, pcr_error, sequencing_error,
                              n_cpu=n_cpu, seeds=seeds)

    elif task == "amplification_rate":
        a = [x/100 for x in range(10, 101)]
        for success_rate in a:
            simulate_multiple(pool_size, final_pool_size, k, pcr_n,
                              success_rate, pcr_error, sequencing_error,
                              n_cpu=n_cpu, seeds=seeds)
            
    elif task == "pool_size":
        a = [i / 100 for i in range(0, 301, 5)]
        a = sorted(set([int(10**i) for i in a]))
        for pool_size in a:
            simulate_multiple(pool_size, final_pool_size, k, pcr_n,
                              success_rate, pcr_error, sequencing_error,
                              n_cpu=n_cpu, seeds=seeds)

    elif task == "final_pool_size":
        a = [i / 100 for i in range(0, 501, 5)]
        a = sorted(set([int(10**i) for i in a]))
        for final_pool_size in a:
            simulate_multiple(pool_size, final_pool_size, k, pcr_n,
                              success_rate, pcr_error, sequencing_error,
                              n_cpu=n_cpu, seeds=seeds)


def simulate_multiple(pool_size, final_pool_size, k, pcr_n,
                      success_rate, pcr_error, sequencing_error,
                      n_cpu=16, seeds=range(10000)):
    '''Instead of accepting seed in simulate(), this function accepts a list of 
    of seed.
    '''
    n_seeds = len(seeds)
    args = zip([pool_size] * n_seeds, [final_pool_size] * n_seeds,
               [k] * n_seeds, [pcr_n] * n_seeds,
               [success_rate] * n_seeds, [pcr_error] * n_seeds,
               [sequencing_error] * n_seeds, seeds)
    with mp.Pool(n_cpu) as p:
        ret = p.starmap(simulate, args)
        for i in ret:
            print(i, end="")

    # for i in seeds:
    #     simulate(pool_size, final_pool_size, k, pcr_n,
    #              success_rate, pcr_error, sequencing_error, seed=i)
    
    
def simulate(pool_size, final_pool_size, k, pcr_n,
             success_rate, pcr_error, sequencing_error, seed):
    '''Simulates one locus'''
    random.seed(a=seed)
    print2("Total number of molecules to start with %d" % pool_size)
    pool = []
    for j in range(pool_size):
        tmp = ''.join(random.choice(('A', 'C', 'G', 'T')) for _ in range(k))
        pool.append(UMIRead("foo", tmp))
    for i in range(pcr_n):
        new_pool = []
        for p in pool:
            new_pool.extend(p.amplify(success_rate, pcr_error))
        print2("Done PCR cycle %d." % (i+1))
        print2("Current number of molecules %d." % (len(new_pool)))
        # This prevents the pool from explosion
        if len(new_pool) <= MAX_POOL_SIZE:
            pool = new_pool
        else:
            pool = random.sample(new_pool, MAX_POOL_SIZE)
            
    print2("Number of reads after PCR: %d" % len(pool))
    if len(pool) > final_pool_size:
        pool2 = random.sample(pool, final_pool_size)
    else:
        pool2 = pool
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
    # print("%s\t%s\t%s" % ("initial_true_n_fragments", "final_n_reads", 
    #                       "final_n_reads_unique", "final_n_reads_1err", 
    #                       "final_n_reads_2err"))

    ########################################################################
    # Allow 0 errors in UMIs. "Unique" method: len(umis)
    ########################################################################
    # dict from observed umi to umi count
    umis = {}
    # dict from true umi to umi count
    umis_true = {}
    for i in final_pool:
        umi = i.seq
        if umi in umis:
            umis[umi] += 1
        else:
            umis[umi] = 1

        umi_true = i.species
        if umi_true in umis_true:
            umis_true[umi_true] += 1
        else:
            umis_true[umi_true] = 1
            
    ########################################################################
    # Allow errors in UMIs. Graph (network) method
    G = umi_graph.UmiGraph(umis, max_ed=1)
    repr_umis = G.get_repr_umi()
    ########################################################################

    ret = "\t".join(("#starting_molecule", "n_reads", "n_umi_true",
                     "n_umi_unique", "n_umi_1err",
                     "pcr_cycle", "umi_length",
                     "initial_pool_size", "final_pool_size",
                     "min_amplification_rate", "pcr_error",
                     "sequencing_error", "seed"))
    my = (pool_size, len(final_pool), len(umis_true),
          len(umis), len(repr_umis),
          pcr_n, k,
          pool_size, final_pool_size,
          success_rate, pcr_error,
          sequencing_error, seed)
    my_str = [str(i) for i in my]
    
    ret2 = "\t".join(my_str)
    return ret + "\n" + ret2 + "\n"

if __name__ == "__main__":
    # test2()
    main()
