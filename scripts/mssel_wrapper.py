import pandas as pd
import sys
import os
import argparse
import subprocess
from math import ceil
import pdb

def writeTraj(df, filename, Ne, npop, ntraj=1):
	
	with open(filename,'w') as traj_File:
		traj_File.write(f"ntraj: {ntraj}\n")
		traj_File.write(f"npop: {npop}\n")

		# print("Running: ploidy = 2, s = " + str(s) + ", N = " + str(N) + ", Dominance = " + str(dom.strip()) + ", Rep = " + str(rep))
		# Sorting was done wrong in original simulation.  Need to sort frequency but leave generations in place
		df1=df["freq"]
		df1=df1.sort_values(ascending=False)
		df1 = pd.concat([df["gen"].reset_index(drop=True), df1.reset_index(drop=True)], axis=1)
		df1['gen'] = df1['gen'].astype(float)
		df1.gen /= float(Ne)
		maxGen = df1['gen'].max()
		df1=df1.loc[(df1.freq.diff() < 0.15),] #any row that differs from previous row by more than 0.15 is likely an error...noticed several odd data points when visualizing data in R.
		if len(df1.loc[(df1.freq.diff() > 0.15),].index) > 3: # Should not find more than a couple artifactual data points
			df1['dif']= df1.freq.diff()
			print(df1)
		traj_File.write(f"n: {len(df1.index) + 1}\t{filename}\n")
		traj_File.write(df1.to_string(index=False,header=None) + "\n")
		traj_File.write(f"{maxGen + 0.000000001}  0  0")
	return(maxGen)

def makeFileNames(outdir, p, s, Ne, dom, rep):
	traj_file = f"{outdir}mssel_input/traj_p{p}_s{s}_N{N}_{dom.strip()}_rep{rep}.txt"
	outfile = f"{outdir}mssel_output/mssel_out_p{p}_s{s}_N{N}_{dom.strip()}_rep{rep}.txt"
	sumout = f"{outdir}mssel_output/mssel_out_p{p}_s{s}_N{N}_{dom.strip()}_rep{rep}.smry.txt"
	return([traj_file, outfile, sumout])

def make2popCmds(mssel_path, n_anc, n_der, p, num_reps, traj_file, selection_spot, rho, n_sites, theta, fuseTime, ms_out, Rscript_path, Rout, num_wind, slide_rate):
	msCmd = f"{mssel_path} { (n_anc + n_der) * p } {num_reps} {n_anc * p} {n_der * p} {traj_file} {selection_spot} -r {rho * p} {n_sites} -t {theta * p} -I 2 0 {n_der * p} {n_anc * p} 0 -ej {fuseTime} 2 1 > {ms_out}"
	rCmd = f"{Rscript_path} {ms_out} {Rout.replace('.txt', '.tmp.txt')} {n_sites} {num_wind} {slide_rate} 2 1 {n_anc * p} {n_der * p}"
	awkCmd = f"awk '{{printf \" {p} {rep} {s} {N} {dom} {n_sites} {theta * p} {rho * p} {Ne} %s {n_anc} {n_der} \\n\", $0}}\' {Rout.replace('.txt', '.tmp.txt')} > {Rout}"
	return([msCmd, rCmd, awkCmd])

# Specify arguments to be read from the command line
parser = argparse.ArgumentParser(description='This script will take a trajectory file produced from Ploidy_Forward_Sim.py, run mssel, and parse the output')
parser.add_argument('-t', type=str, metavar='trajectory_file', required=True, help='Full path to results from Ploidy_Forward_Sim.py')
parser.add_argument('-o', type=str, metavar='output_directory', required=True, help='Full path to output directory.')
parser.add_argument('-ms', type=str, metavar='mssel_path', default="/Users/pmonnahan/Documents/Research/code/dmc/mssel_modified/mssel", required=False, help='Full path to mssel executable')
parser.add_argument('-X', type=str, metavar='summary_script', default="/Users/pmonnahan/Documents/Research/PloidySim/scripts/ms_summarize.R", required=False, help='Path to ms_summarize.R')
parser.add_argument('-reps', type=int, metavar='number_of_reps', default="1", required=False, help='Number of mssel reps for a single trajectory')
parser.add_argument('-mu', type=float, metavar='mutation_rate', default="0.000000015", required=False, help='specified as decimal.')
parser.add_argument('-r', type=float, metavar='recombination_rate', default="0.000000015", required=False, help='specified as decimal.')
parser.add_argument('-L', type=int, metavar='sequence_length', default="200000", required=False, help='Total number of sites')
parser.add_argument('-T', type=float, metavar='Theta', default=-9.0, required=False, help='If provided, this will override calculation of theta from pop size under simulation')
parser.add_argument('-R', type=float, metavar='Rho', default=-9.0, required=False, help='If provided, this will override calculation of rho from pop size under simulation')
parser.add_argument('-na', type=int, metavar='n_ancestral', default=10, required=False, help='number of individuals possessing the ancestral (i.e non-sweep) allele')
parser.add_argument('-nd', type=int, metavar='n_derived', default=10, required=False, help='number of individuals possessing the derived (i.e. sweep) allele')
parser.add_argument('-w', type=str, metavar='number_of_windows', default="200", required=False, help='number of windows to be used in ms_summarize.R')
parser.add_argument('-W', type=str, metavar='window_slide_rate', default="0.5", required=False, help='Proportion of winSize that we want to slide the window')
parser.add_argument('-Ne', type=int, metavar='effective_pop_size', default=-9, required=False, help='use a single value of Ne instead of using the N under which sweeps were simulated')
args = parser.parse_args()

args.o += "/"
if os.path.exists(args.o) is False: os.mkdir(args.o)
if os.path.exists(args.o + "mssel_input") is False: os.mkdir(args.o + "mssel_input")
if os.path.exists(args.o + "mssel_output") is False: os.mkdir(args.o + "mssel_output")

outdir = args.o
mssel_path = args.ms
Rscript_path = args.X
mu = args.mu
r = args.r
num_reps = args.reps
num_wind = args.w
slide_rate = args.W
n_anc = args.na # number of alleles sampled from ancestral population
n_der = args.nd # number of alleles sampled from selected population
n_sites = args.L #number of sites
selection_spot = ceil(n_sites * 0.5)
segsites=1000
npop = 2
GPSS = 100 #Generations prior to start of selection

#Is Pandas the way to go on this? Takes a long time to read data into memory
TRAJ = pd.read_csv(args.t, header=0, dtype={'rep': object, 'ploidy': float, 'start': float, 'N': float, 's': float, 'dom': object, 'gen': float, 'Ngen': float, 'freq': float})

# Loop over
jobList=[] 
for idx, df in TRAJ.groupby(['ploidy', 's', 'N', 'dom', 'rep']):

	p = int(df['ploidy'].iloc[0])
	if p == 2:
		alt_fact = 2
	else:
		alt_fact = 0.5

	s = float(df['s'].iloc[0])
	N = int(df['N'].iloc[0])
	dom = str(df['dom'].iloc[0])
	rep = str(df['rep'].iloc[0])

	if args.Ne == -9: Ne = N
	else: Ne = args.Ne
	if args.T != -9.0:
		theta = args.T
		# n_sites = ceil(theta / (2 * mu * p * Ne))
		n_sites = ceil(theta / (2 * mu * Ne))
		selection_spot = ceil(n_sites * 0.5)
	else: theta = ceil(2 * Ne * mu * n_sites)
	if args.R != -9.0: rho = args.R
	else: rho = ceil(2 * Ne * r * n_sites)

	for dd in range(0, num_reps):
		rep1 = f"{rep}-{dd}"
		files = makeFileNames(outdir, p, s, Ne, dom, rep1)
		maxGen = writeTraj(df, files[0], Ne * p, npop)
		fuseTime = maxGen + (GPSS / (Ne * p))

		#Delete intermediate files
		cmds = make2popCmds(mssel_path, n_anc, n_der, p, "1", files[0], selection_spot, rho, n_sites, theta, fuseTime, files[1], Rscript_path, files[2], num_wind, slide_rate)
		
		print(";".join(cmds))

	p *= alt_fact

	for dd in range(0, num_reps):
		rep2 = f"{rep}-{dd}"
		files = makeFileNames(outdir, p, s, Ne, dom, f"{rep2}Alt")
		maxGen = writeTraj(df, files[0], Ne * p, npop)
		fuseTime = maxGen + (GPSS / (Ne * p))

		alt_cmds = make2popCmds(mssel_path, n_anc, n_der, p, "1", files[0], selection_spot, rho, n_sites, theta, fuseTime, files[1], Rscript_path, files[2], num_wind, slide_rate)

		print(";".join(alt_cmds))