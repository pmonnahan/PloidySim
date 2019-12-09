import pdb
import random
import argparse
import numpy
import copy

def initialize(size, ploidy, freq): 
    pop = [[] for x in range(0, size)] #Each individual is a list of alleles (0 or 1) and the population is a list of lists.
    for i, ind in enumerate(pop):
        for k in range(0,ploidy):
            rx = random.random()
            if rx < freq: # Draw alleles for each individual based on starting frequency
                pop[i].append(1)
            else:
                pop[i].append(0)
    return(pop)

def drift(pop, numGens):
    popSize = len(pop)
    p = calcFreq(pop) 
    F = []
    for i in range(0,numGens): #Simulate drift for specified number of generations
        if p > 0 and p < 1:
            newpop = [] #Store individuals for next generation's population as they are created
            for j in range(0, popSize):
                p1 = copy.deepcopy(random.choice(pop)) #Parent 1 is a random individual drawn from population
                p2 = copy.deepcopy(random.choice(pop)) #Parent 2

                g1 = [p1.pop() for x in range(0, int(len(p1)/2))] #Gamete1. Randomly sample alleles without replacement from Parent1
                g2 = [p2.pop() for x in range(0, int(len(p2)/2))] #Gamete2

                o1 = g1 + g2 #New individual is sum of gamete lists
                newpop.append(o1) #Add new individual to next generations' population
            p = calcFreq(newpop)
            F.append(p)

        # print(i, p)
        pop = newpop
        
    return(pop, F)

def calcFreq(pop):
    M = 0
    T = 0
    for ind in pop:
        M += sum(ind) #Tracks number of mutant (1) alleles in the population
        T += len(ind) #Tracks total alleles
    if T != 0: p = M / T
    else: p = 0
    return(p)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-N', type=int, metavar='Population_Size', required=False, default=1000)
    parser.add_argument('-p', type=float, metavar='start_freq', required=False, default=0.5)
    parser.add_argument('-k', type=int, metavar='ploidy', required=False, default = 2)
    parser.add_argument('-t', type=int, metavar='numGens', required=False, default=100)
    parser.add_argument('-r', type=int, metavar='replicates', required=False, default = 1000)

    # parser.add_argument('-r', action='store_true', help='restrictTo_mainChroms')

    args = parser.parse_args()

    start_freq = args.p
    N = args.N
    t = args.t 
    k = args.k #ploidy
    reps = args.r

    for k in [2, 4, 8]:
        freq_diffs = []
        varps = [] 
        for rep in range(0, reps):
            pop = initialize(N, k, start_freq)
            start_freq = calcFreq(pop)
            if start_freq > 0 and start_freq < 1:
                newpop, traj = drift(pop, t)
                end_freq = calcFreq(newpop)
                freq_diffs.append(abs(start_freq - end_freq))

                dp = [x - y for x, y in zip(traj, traj[1:])] #Calculate Delta p in each generation based on drift trajectory
                if numpy.var(dp) != "nan":
                    print(f"{N} {k} {start_freq} {rep} {t} {numpy.var(dp)}") #Drift rate is measured as variance in allele frequency change in each generation
                    varps.append(numpy.var(dp))


        # for gen, freq in enumerate(dp):
        #     print(f"{N} {k} {start_freq} {rep} {t} {gen} {freq}")


    # print(numpy.mean(freq_diffs))
    # print(numpy.nanmean(varps))     