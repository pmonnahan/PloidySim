import pandas as pd
import sys
import os
import argparse
import subprocess
from math import ceil

# Specify arguments to be read from the command line
parser = argparse.ArgumentParser(description='This script will take a trajectory file produced from Ploidy_Forward_Sim.py, run mssel, and parse the output')
parser.add_argument('-t', type=str, metavar='trajectory_file', required=True, help='Full path to results from Ploidy_Forward_Sim.py')
parser.add_argument('-o', type=str, metavar='output_directory', required=True, help='Full path to output directory.')
parser.add_argument('-ms', type=str, metavar='mssel_path', default="/Users/pmonnahan/Documents/Research/code/dmc/mssel_modified/mssel", required=False, help='Full path to mssel executable')
parser.add_argument('-X', type=str, metavar='summary_script', default="/Users/pmonnahan/Documents/Research/PloidySim/scripts/ms_summarize.R", required=False, help='Path to ms_summarize.R')
parser.add_argument('-reps', type=str, metavar='number_of_reps', default="1", required=False, help='Number of mssel reps for a single trajectory')
parser.add_argument('-mu', type=float, metavar='mutation_rate', default="0.000000015", required=False, help='specified as decimal.')
parser.add_argument('-r', type=float, metavar='recombination_rate', default="0.000000015", required=False, help='specified as decimal.')
parser.add_argument('-L', type=int, metavar='sequence_length', default="200000", required=False, help='Total number of sites')
parser.add_argument('-T', type=float, metavar='Theta', default=-9.0, required=False, help='If provided, this will override calculation of theta from pop size under simulation')
parser.add_argument('-R', type=float, metavar='Rho', default=-9.0, required=False, help='If provided, this will override calculation of rho from pop size under simulation')
parser.add_argument('-na', type=int, metavar='n_ancestral', default=10, required=False, help='number of individuals possessing the ancestral (i.e non-sweep) allele')
parser.add_argument('-nd', type=int, metavar='n_derived', default=10, required=False, help='number of individuals possessing the derived (i.e. sweep) allele')
parser.add_argument('-w', type=str, metavar='window_size', default="50", required=False, help='window size in number of snps for ms_summarize.R')
parser.add_argument('-W', type=str, metavar='window_slide_rate', default="0.5", required=False, help='Proportion of winSize that we want to slide the window')
args = parser.parse_args()

args.o += "/"
if os.path.exists(args.o) is False: os.mkdir(args.o)
if os.path.exists(args.o + "mssel_input") is False: os.mkdir(args.o + "mssel_input")
if os.path.exists(args.o + "mssel_output") is False: os.mkdir(args.o + "mssel_output")

mu = args.mu
r = args.r
num_reps = args.reps
n_anc = args.na # number of alleles sampled from ancestral population
n_der = args.nd # number of alleles sampled from selected population
n_sites = args.L #number of sites
selection_spot = ceil(n_sites * 0.5)
segsites=1000

#Is Pandas the way to go on this? Takes a long time to read data into memory
TRAJ = pd.read_csv(args.t, header=0, dtype={'rep': object, 'ploidy': float, 'start': float, 'N': float, 's': float, 'dom': object, 'gen': float, 'Ngen': float, 'freq': float})

# Loop over
jobList=[] 
for idx, df in TRAJ.groupby(['ploidy', 's', 'N', 'dom', 'rep']):
	p = int(df['ploidy'].iloc[0])
	s = float(df['s'].iloc[0])
	N = int(df['N'].iloc[0])
	dom = str(df['dom'].iloc[0])
	rep = str(df['rep'].iloc[0])

	traj_file = f"{args.o}mssel_input/traj_p{p}_s{s}_N{N}_{dom.strip()}_rep{rep}.txt"
	traj_File=open(traj_file,'w')
	traj_File.write("ntraj: 1\n")
	traj_File.write("npop: 1\n")
	
	outfile = f"{args.o}mssel_output/mssel_out_p{p}_s{s}_N{N}_{dom.strip()}_rep{rep}.txt"
	sumout = f"{args.o}mssel_output/mssel_out_p{p}_s{s}_N{N}_{dom.strip()}_rep{rep}.smry.txt"
	print("Running: ploidy = 2, s = " + str(s) + ", N = " + str(N) + ", Dominance = " + str(dom.strip()) + ", Rep = " + str(rep))
	df1=df.loc[df.rep == rep,("Ngen","freq")]
	df1['Ngen'] = df1['Ngen'].astype(float)
	df1=df1.sort_values(by=["Ngen"])
	df1=df1.loc[(df1.freq.diff() < 0.15),] #any row that differs from previous row by more than 0.15 is likely an error...noticed several odd data points when visualizing data in R.
	if len(df1.loc[(df1.freq.diff() > 0.15),].index) > 3: # Should not find more than a couple artifactual data points
		df1['dif']= df1.freq.diff()
		print(df1)
	traj_File.write(f"n: {len(df1.index)}\trep{rep}\n")
	traj_File.write(df1.to_string(index=False,header=None) + "\n")

	if args.T != -9.0:
		theta = args.T
		n_sites = ceil(theta / (2 * mu * p * N))
		selection_spot = ceil(n_sites * 0.5)
	else:
		theta = 2 * p * N * mu * n_sites
	if args.R != -9.0:
		rho = args.R
	else:
		rho = 2 * p * N * r * n_sites

	msCmd = f"{args.ms} { (n_anc + n_der) * p } {num_reps} {n_anc * p} {n_der * p} {traj_file} {selection_spot} -r {rho} {n_sites} -t {theta}  > {outfile}; "
	rCmd = f"{args.X} {outfile} {sumout.replace('.txt', '.tmp.txt')} {n_sites} {args.w} {args.W} 2 1 {n_anc} {n_der}; "
	awkCmd = f"awk '{{printf \" {p} {rep} {s} {N} {dom} %s {n_anc} {n_der} \\n\", $0}}\' {sumout.replace('.txt', '.tmp.txt')} > {sumout}"
	
	print(msCmd + rCmd + awkCmd)
				# pp = subprocess.Popen(cmd,shell = True)

# for s in tetTraj.s.unique():
# 	for N in tetTraj.N.unique():
# 		for dom in tetTraj.dom.unique():
# 			df=tetTraj.loc[(tetTraj.s == s) & (tetTraj.N == N) & (tetTraj.dom == dom),]
# 			for rep in df.rep.unique():
# 				traj_file = args.o + "mssel_input/traj_p4_s" + str(s) + "_N" + str(N) + "_" + dom.strip() + "_rep" + str(rep) + ".txt "
# 				traj_File=open(traj_file,'w')
# 				traj_File.write("ntraj: " + str(len(traj.rep.unique()) - 1) + "\n")
# 				traj_File.write("npop: 1\n")
				
# 				outfile = args.o + "mssel_output/mssel_out_p4_s" + str(s) + "_N" + str(N) + "_" + dom.strip() + "_rep" + str(rep) + ".txt"
# 				# print("Running: ploidy = 2, s = " + str(s) + ", N = " + str(N) + ", Dominance = " + str(dom.strip()) + ", Rep = " + str(rep))
# 				df1=df.loc[df.rep == rep,("Ngen","freq")]
# 				df1['Ngen'] = df1['Ngen'].astype(float)
# 				df1=df1.sort_values(by=["Ngen"])
# 				df1=df1.loc[(df1.freq.diff() < 0.15),] #any row that differs from previous row by more than 0.15 is likely an error...noticed several odd data points when visualizing data in R.
# 				if len(df1.loc[(df1.freq.diff() > 0.15),].index) > 3: # Should not find more than a couple artifactual data points
# 					df1['dif']= df1.freq.diff()
# 					print(df1)
# 				traj_File.write("n: " + str(len(df1.index)) + "\trep" + str(rep) + "\n")
# 				traj_File.write(df1.to_string(index=False,header=None) + "\n")

# 				theta = 8 * N * mu * n_sites
# 				rho = 8 * N * r * n_sites

# 				# To end of cmd, add && R analysis to summarise the data
# 				cmd = args.ms + ' ' + str((n_anc * 4) + (n_der * 4)) + ' ' + str(num_reps) + ' ' + str(n_anc * 4) + ' ' + str(n_der * 4) + ' ' + traj_file + str(selection_spot) + ' -r ' + str(rho) + ' ' + str(n_sites) + ' -s ' + str(segsites) + " > " + outfile
# 				print(cmd)
# 				# pp = subprocess.Popen(cmd,shell = True)
