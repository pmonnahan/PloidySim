import pandas as pd
import sys
import os
import argparse
import subprocess

# Specify arguments to be read from the command line
parser = argparse.ArgumentParser(description='This script will take a trajectory file produced from Ploidy_Forward_Sim.py, run mssel, and parse the output')
parser.add_argument('-t', type=str, metavar='trajectory_file', required=True, help='Full path to results from Ploidy_Forward_Sim.py')
parser.add_argument('-o', type=str, metavar='output_directory', required=True, help='Full path to output directory.')
parser.add_argument('-ms', type=str, metavar='mssel_path', default="/Users/pmonnahan/Documents/Research/code/dmc/mssel_modified/mssel", required=False, help='Full path to mssel executable')
args = parser.parse_args()

args.o += "/"
if os.path.exists(args.o) is False: os.mkdir(args.o)
if os.path.exists(args.o + "mssel_input") is False: os.mkdir(args.o + "mssel_input")
if os.path.exists(args.o + "mssel_output") is False: os.mkdir(args.o + "mssel_output")

mu = 10**(-8)
r = 10**(-8)
num_reps = 10
n_anc = 0 # number of alleles sampled from ancestral population
n_der = 10 # number of alleles sampled from selected population
n_sites = 200000 #number of sites
selection_spot = 25000
segsites=1000

#Is Pandas the way to go on this? Takes a long time to read data into memory
TRAJ = pd.read_csv(args.t, header=0)

samp_sizes = args.n.split(",")
# Loop over
jobList=[] 
for idx, df in TRAJ.groupby(['ploidy', 's', 'N', 'dom', 'rep']):
	p = int(df['ploidy'].iloc[0])
	s = str(df['s'].iloc[0])
	N = int(df['N'].iloc[0])
	dom = str(df['dom'].iloc[0])
	rep = str(df['rep'].iloc[0])

	traj_file = f"{args.o}mssel_input/traj_p{p}_s{s}_N{N}_{dom.strip()}_rep{rep}.txt"
	traj_File=open(traj_file,'w')
	traj_File.write("ntraj: 1\n")
	traj_File.write("npop: 1\n")
	
	outfile = f"{args.o}mssel_output/mssel_out_p{p}_s{s}_N{N}_{dom.strip()}_rep{rep}.txt"
	# print("Running: ploidy = 2, s = " + str(s) + ", N = " + str(N) + ", Dominance = " + str(dom.strip()) + ", Rep = " + str(rep))
	df1=df.loc[df.rep == rep,("Ngen","freq")]
	df1['Ngen'] = df1['Ngen'].astype(float)
	df1=df1.sort_values(by=["Ngen"])
	df1=df1.loc[(df1.freq.diff() < 0.15),] #any row that differs from previous row by more than 0.15 is likely an error...noticed several odd data points when visualizing data in R.
	if len(df1.loc[(df1.freq.diff() > 0.15),].index) > 3: # Should not find more than a couple artifactual data points
		df1['dif']= df1.freq.diff()
		print(df1)
	traj_File.write(f"n: {len(df1.index)}\trep{rep}\n")
	traj_File.write(df1.to_string(index=False,header=None) + "\n")

	theta = 2 * p * N * mu * n_sites
	rho = 2 * p * N * r * n_sites

	cmd = f"{args.ms} {(n_anc + n_der) * p} + ' ' + str(num_reps) + ' ' + str(n_anc * 2) + ' ' + str(n_der * 2) + ' ' + traj_file + str(selection_spot) + ' -r ' + str(rho) + ' ' + str(n_sites) + ' -s ' + str(segsites) + " > " + outfile
	print(cmd)
				# pp = subprocess.Popen(cmd,shell = True)

for s in tetTraj.s.unique():
	for N in tetTraj.N.unique():
		for dom in tetTraj.dom.unique():
			df=tetTraj.loc[(tetTraj.s == s) & (tetTraj.N == N) & (tetTraj.dom == dom),]
			for rep in df.rep.unique():
				traj_file = args.o + "mssel_input/traj_p4_s" + str(s) + "_N" + str(N) + "_" + dom.strip() + "_rep" + str(rep) + ".txt "
				traj_File=open(traj_file,'w')
				traj_File.write("ntraj: " + str(len(traj.rep.unique()) - 1) + "\n")
				traj_File.write("npop: 1\n")
				
				outfile = args.o + "mssel_output/mssel_out_p4_s" + str(s) + "_N" + str(N) + "_" + dom.strip() + "_rep" + str(rep) + ".txt"
				# print("Running: ploidy = 2, s = " + str(s) + ", N = " + str(N) + ", Dominance = " + str(dom.strip()) + ", Rep = " + str(rep))
				df1=df.loc[df.rep == rep,("Ngen","freq")]
				df1['Ngen'] = df1['Ngen'].astype(float)
				df1=df1.sort_values(by=["Ngen"])
				df1=df1.loc[(df1.freq.diff() < 0.15),] #any row that differs from previous row by more than 0.15 is likely an error...noticed several odd data points when visualizing data in R.
				if len(df1.loc[(df1.freq.diff() > 0.15),].index) > 3: # Should not find more than a couple artifactual data points
					df1['dif']= df1.freq.diff()
					print(df1)
				traj_File.write("n: " + str(len(df1.index)) + "\trep" + str(rep) + "\n")
				traj_File.write(df1.to_string(index=False,header=None) + "\n")

				theta = 8 * N * mu * n_sites
				rho = 8 * N * r * n_sites

				# To end of cmd, add && R analysis to summarise the data
				cmd = args.ms + ' ' + str((n_anc * 4) + (n_der * 4)) + ' ' + str(num_reps) + ' ' + str(n_anc * 4) + ' ' + str(n_der * 4) + ' ' + traj_file + str(selection_spot) + ' -r ' + str(rho) + ' ' + str(n_sites) + ' -s ' + str(segsites) + " > " + outfile
				print(cmd)
				# pp = subprocess.Popen(cmd,shell = True)
