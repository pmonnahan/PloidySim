import collections
import numpy as np
import argparse


def simulate(ploidy, starting_freq, pop_size, s, dom):

    pop = np.random.binomial(ploidy, starting_freq, pop_size)

    # Makes a dictionary with keys equal to allele count and value equal to number of individuals with this allele count
    counter = collections.Counter(pop)

    geno_freq = [0.0 for x in range(0, ploidy + 1)]
    # can we replace g_prime with p_prime?  
    counts = list(counter.values())
    indexes = list(counter.keys())
    for x in indexes:
        geno_freq[x] = counts[x] / pop_size

    af = sum(pop) / (ploidy * pop_size)
    gen = 1
    Ngen = gen / (2 * ploidy * pop_size)

    fitness = [0 for x in range(0, ploidy + 1)]
    fitness[ploidy] = 1
    fitness[0] = 1 - s

    for f in range(1, ploidy):
        if dom == "additive":
            h = s / ploidy
            fitness[f] = (1 - s) + (f * h) #Is this valid?
        elif dom == "recessive":
            fitness[f] = (1 - s)    
        elif dom == "dominant":
            fitness[f] = 1

    af_list = [af]
    Ngen_list = [Ngen]
    while af < 0.99:
        counter = collections.Counter(pop)

        geno_freq = [0.0 for x in range(0, ploidy + 1)]
        # can we replace g_prime with p_prime?  
        counts = list(counter.values())
        indexes = list(counter.keys())
        # print(counts,indexes)
        # print(geno_freq)

        Counts = [0 for i in range(0, ploidy + 1)]
        for j, i in enumerate(indexes):
            Counts[i] = counts[j]
        geno_freq = [x / pop_size for x in Counts]

        tot_fitness = sum([a * b for a, b in zip(geno_freq, fitness)])
        g_prime = [a * b / tot_fitness for a, b in zip(geno_freq, fitness)] 
        # print(af)
        # print(geno_freq)
        # print(g_prime)
        if af == 0.0:
            break
        # Make individuals for next generation
        new_pop = []
        for n in range(0, pop_size):

            rx1 = np.random.random()
            rx2 = np.random.random()
            thresh = 0
            p1 = -9
            p2 = -9
            # Draw parents
            for i, g in enumerate(g_prime):
                if rx1 > thresh and rx1 <= thresh + g:
                    p1 = i
                if rx2 > thresh and rx2 <= thresh + g:
                    p2 = i
                thresh += g

            assert (p1 != -9 and p2 != -9)

            # Make gametes from selected parents and merge to create new individual
            new_ind = 0
            for j in range(0, int(ploidy / 2)):
                af1 = p1 / (ploidy - j)
                af2 = p2 / (ploidy - j)
                rx3 = np.random.random()
                rx4 = np.random.random()
                if rx3 < af1:
                    new_ind += 1
                    p1 -= 1
                if rx4 < af2:
                    new_ind += 1
                    p2 -= 1
            new_pop.append(new_ind)
        pop = new_pop
        gen += 1
        af = sum(pop) / (ploidy * pop_size)
        Ngen = gen / (2 * ploidy * pop_size)
        af_list.append(af)
        Ngen_list.append(Ngen)

    return af_list, Ngen_list


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Forward simulation of selection on beneficial allele.')
    parser.add_argument('-r', type=int, metavar='replicates', required=True, help='number of allele frequency trajectories to report')
    parser.add_argument('-p', type=int, metavar='ploidy', required=True, help='code ploidy as integer 2 for diploids 4 for tetraploids')
    parser.add_argument('-f', type=float, metavar='starting_frequency', required=True, help='value from 0 to 1')
    parser.add_argument('-N', type=int, metavar='Population_Size', required=True, help='Number of individuals in the population.')
    parser.add_argument('-s', type=float, metavar='selection_coefficient', required=True, help='wild-type homozygote has fitness 1-s.')
    parser.add_argument('-d', type=str, metavar='dominance_coefficient', required=True, help='specified for beneficial allele: dominant, recessive, or additive')
    
    args = parser.parse_args()

    for j in range(0, args.r):
        fixed = False
        while fixed is False:
            af, Ngen = simulate(args.p, args.f, args.N, args.s, args.d)
            if af[-1] == 0.0:
                fixed = False
            elif af[-1] >= 0.99:
                fixed = True
            else:
                print("Neither allele fixed")
        for i, k in enumerate(af):
            print(j, "\t", args.p, "\t", args.f, "\t", args.N, "\t", args.s, "\t", args.d, "\t", i, "\t", Ngen[i], "\t", k)

