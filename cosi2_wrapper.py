import pandas as pd
import sys
import os
import argparse
import subprocess

# Specify arguments to be read from the command line
parser = argparse.ArgumentParser(description='This script will take a trajectory file produced from Ploidy_Forward_Sim.py, run mssel, and parse the output')
parser.add_argument('-t', type=str, metavar='trajectory_file', required=True, help='Full path to results from Ploidy_Forward_Sim.py')
parser.add_argument('-o', type=str, metavar='output_directory', required=True, help='Full path to output directory.')
parser.add_argument('-mu', type=str, metavar='mutation_rate', default="0.000000015", required=False, help='specified as decimal.')
parser.add_argument('-r', type=str, metavar='recombination_file', default="/Users/pmonnahan/Documents/Research/PloidySim/cosi2_input/model.test", required=False, help='Recombination map file')
parser.add_argument('-reps', type=str, metavar='number_of_reps', default="1", required=False, help='Number of cosi2 reps for a single trajectory')
parser.add_argument('-n', type=str, metavar='sample_sizes', default="10,10", required=False, help='sample size of each population specified as comma separated list.')
parser.add_argument('-N', type=str, metavar='population_sizes', default="-9", required=False, help='size of each population')
parser.add_argument('-selPop', type=int, metavar='selection_population', default="1", required=False, help='Population in which the sweep will occur; 0-based.')
parser.add_argument('-selStop', type=str, metavar='selection_stop', default="0", required=False, help='number of generations before present at which selection ended')
parser.add_argument('-selSpot', type=str, metavar='selection_spot', default="0.5", required=False, help='Location along sequence of sweeping mutation.  Specified on [0,1]')
parser.add_argument('-L', type=str, metavar='sequence_length', default="200000", required=False, help='Total number of sites')

args = parser.parse_args()

args.o += "/"
if os.path.exists(args.o) is False: os.mkdir(args.o)
if os.path.exists(args.o + "cosi2_input") is False: os.mkdir(args.o + "cosi2_input")
if os.path.exists(args.o + "cosi2_output") is False: os.mkdir(args.o + "cosi2_output")


#Is Pandas the way to go on this? Takes a long time to read data into memory
traj = pd.read_table(args.t,sep=r"\s*",header=0)

dipTraj = traj.loc[traj['ploidy'] == 2,]
tetTraj = traj.loc[traj['ploidy'] == 4,]

# Loop over
jobList=[] 
for s in dipTraj.s.unique():
	for N in dipTraj.N.unique():
		for dom in dipTraj.dom.unique():
			df=dipTraj.loc[(dipTraj.s == s) & (dipTraj.N == N) & (dipTraj.dom == dom),]
			for rep in df.rep.unique():

				#These files will hold the sweep trajectory for input into cosi2
				traj_file = args.o + "cosi2_input/traj_p2_s" + str(s) + "_N" + str(N) + "_" + dom.strip() + "_rep" + str(rep) + ".txt "				
				
				# print("Running: ploidy = 2, s = " + str(s) + ", N = " + str(N) + ", Dominance = " + str(dom.strip()) + ", Rep = " + str(rep))
				df1=df.loc[df.rep == rep,("gen","freq")]
				df1['gen'] = df1['gen'].astype(int) + int(args.selStop)
				df1=df1.sort_values(by=["gen"])
				df1=df1.loc[(df1.freq.diff() < 0.15),] #any row that differs from previous row by more than 0.15 is likely an error...noticed several odd data points when visualizing data in R.
				if len(df1.loc[(df1.freq.diff() > 0.15),].index) > 3: # Should not find more than a couple artifactual data points
					df1['dif']= df1.freq.diff()
					print(df1)
				
				df1.to_csv(traj_file, index=False,header=None, sep=" ")

				## VERIFY EVERYTHING...THETA IS 4 N IN DIPS AND 8 N IN TETS?
				theta = 4 * N * float(args.mu) * float(args.L)

				# Write parameter file
				param_file = args.o + "cosi2_input/param_p2_s" + str(s) + "_N" + str(N) + "_" + dom.strip() + "_rep" + str(rep) + ".txt "
				param_File = open(param_file,'w')
				param_File.write("# p = 2; s = " + str(s) + "; N = " + str(N) + "\n")
				param_File.write("length " + args.L + "\n\n")
				param_File.write("mutation_rate " + args.mu + "\n\n")
				param_File.write("recomb_file " + args.r + "\n\n")

				#Define populations
				for i, samp in enumerate(args.n.split(",")):
					if i == args.selPop:
						param_File.write("pop_define " + str(i + 1) + " sweep\n")
					else:
						param_File.write("pop_define " + str(i + 1) + " standard\n")
					if args.N != "-9":
						pop_sizes = args.N.split(",")
						assert len(pop_sizes) >= args.selPop + 1
						param_File.write("pop_size " + str(i + 1) + " " + pop_sizes[i])
					else:
						param_File.write("pop_size " + str(i + 1) + " " + str(N))
					param_File.write("sample_size " + str(i + 1) + " " + samp)

				#Define population events

				outfile = args.o + "cosi2_output/cosi2_out_p2_s" + str(s) + "_N" + str(N) + "_" + dom.strip() + "_rep" + str(rep) + ".txt"
				cmd = "env COSI_NEWSIM=1 COSI_LOAD_TRAJ=" + traj_file + " coalescent -p" + param_file + " -m > " + outfile
				
				# pp = subprocess.Popen(cmd,shell = True)




