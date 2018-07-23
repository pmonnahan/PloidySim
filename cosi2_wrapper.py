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
parser.add_argument('-L', type=str, metavar='sequence_length', default="200000", required=False, help='Total number of sites')
parser.add_argument('-events', type=str, metavar='population_events', required=True, help='comma-separated list of population events (e.g. split "pop1/pop2 split" 2 1 B100,sweep "sweep1" 2 256 0.1 0.5 0.95). When specifying time values, include B, D, or A in front of the value (e.g. B100).  B specifies that the event will take place X generations before selection begins. D (during) specifies that the event happens X generations before selection ends, and A specifies that the event takes place after selection has ended.  For sweeps, the specified time indicates the number of generations that selection ended prior to present time.  Aside from these time values, population event specification follows the cosi2 manual at https://software.broadinstitute.org/mpg/cosi2/cosidoc.html')
parser.add_argument('-suf', type=str, metavar='param_file_suffix', default="", required=False, help='Suffix to add to end of param_file for accounting purposes')
parser.add_argument('-X', type=str, metavar='summary_script', default="/Users/pmonnahan/Documents/Research/PloidySim/ms_summarize.R", required=False, help='Path to ms_summarize.R')
parser.add_argument('-w', type=str, metavar='window_size', default="50", required=False, help='window size in number of snps for ms_summarize.R')
parser.add_argument('-W', type=str, metavar='window_slide_rate', default="0.5", required=False, help='Proportion of winSize that we want to slide the window')

args = parser.parse_args()

args.o += "/"
if os.path.exists(args.o) is False: os.mkdir(args.o)
if os.path.exists(args.o + "cosi2_input") is False: os.mkdir(args.o + "cosi2_input")
if os.path.exists(args.o + "cosi2_output") is False: os.mkdir(args.o + "cosi2_output")


#Is Pandas the way to go on this? Takes a long time to read data into memory
traj = pd.read_table(args.t,sep=r"\s*",header=0)

dipTraj = traj.loc[traj['ploidy'] == 2,]
tetTraj = traj.loc[traj['ploidy'] == 4,]

samp_sizes = args.n.split(",")
num_pops = len(samp_sizes)

# Loop over
jobList=[] 
for s in dipTraj.s.unique():
	for N in dipTraj.N.unique():
		for dom in dipTraj.dom.unique():
			df=dipTraj.loc[(dipTraj.s == s) & (dipTraj.N == N) & (dipTraj.dom == dom),]
			for rep in df.rep.unique():

				# find generation at which selection stopped from the sweep pop_event and the population in which the sweep will occur
				events = args.events.split(",")
				sweeps = [s for s in events if "sweep" in s]
				assert len(sweeps) == 1, "Multiple sweeps specified in events.  Not currently supported."
				selStop = sweeps[0].split()[3]
				selPop = sweeps[0].split()[2]

				#Retrieve and format data for each rep of a set of run parameters
				# print("Running: ploidy = 2, s = " + str(s) + ", N = " + str(N) + ", Dominance = " + str(dom.strip()) + ", Rep = " + str(rep))
				df1=df.loc[df.rep == rep,("gen","freq")]
				selGen = len(df1.index) # Number of generations over which selection will occur
				df1['gen'] = abs(df1['gen'].astype(int) - selGen) + int(selStop)
				df1=df1.sort_values(by=["gen"])
				df1=df1.loc[(df1.freq.diff() < 0.15),] #any row that differs from previous row by more than 0.15 is likely an error...noticed several odd data points when visualizing data in R.
				if len(df1.loc[(df1.freq.diff() > 0.15),].index) > 3: # Should not find more than a couple artifactual data points
					df1['dif']= df1.freq.diff()
					print(df1)

				maxFreq=df1['freq'].max()

				# Write parameter file
				param_file = args.o + "cosi2_input/param_p2_s" + str(s) + "_N" + str(N) + "_" + dom.strip() + "_rep" + str(rep) + "_" + args.suf + ".txt"
				param_File = open(param_file,'w')
				param_File.write("# p = 2; s = " + str(s) + "; N = " + str(N) + "\n")
				param_File.write("length " + args.L + "\n\n")
				param_File.write("mutation_rate " + args.mu + "\n\n")
				param_File.write("recomb_file " + args.r + "\n\n")

				#Define populations
				for i, samp in enumerate(samp_sizes):
					if str(i) == selPop:
						param_File.write("pop_define " + str(i + 1) + " sweep" + str(i + 1) + "\n")
					else:
						param_File.write("pop_define " + str(i + 1) + " nosweep" + str(i + 1) + "\n")
					if args.N != "-9":
						pop_sizes = args.N.split(",")
						assert len(pop_sizes) >= selPop + 1
						param_File.write("pop_size " + str(i + 1) + " " + str(int(pop_sizes[i]) * 2) + "\n")
					
					##### IS THIS VALID? MULTIPLYING POPULATION SIZE BY PLOIDY? ####
					else: 
						param_File.write("pop_size " + str(i + 1) + " " + str(N * 2) + "\n")
					##### NEEDS TO BE CHECKED #####

					param_File.write("sample_size " + str(i + 1) + " " + str(int(samp * 2)) + "\n\n")

				#Define population events
				for event in events:
					event = event.split()
					if event[0] == "sweep":
						event[3] = str(selGen + int(selStop))
						event.append(str(maxFreq))
					elif event[0] == "split":
						print(event[4][1:])
						T = int(event[4][1:])
						if event[4].startswith("B"):
							event[4] = str(T + selGen + int(selStop))
						elif event[4].startswith("D"):	
							event[4] = str(int(selStop) + T)
						elif event[4].startswith("A"):
							assert int(selStop) > T, "Negative (future) time for event: %r" % event[1]
							event[4] = str(int(selStop) - T)
					else:
						print("Event of type '" + event[0] + "' is not currently implemented.")
					param_File.write("pop_event " + " ".join(event) + "\n\n")
				param_File.close()

				#This files will hold the sweep trajectory for input into cosi2
				traj_file = args.o + "cosi2_input/traj_p2_s" + str(s) + "_N" + str(N) + "_" + dom.strip() + "_rep" + str(rep) + "_" + args.suf + ".txt"
				df1.to_csv(traj_file, index=False,header=None, sep=" ")

				#File that will store output of cosi2
				outfile = args.o + "cosi2_output/cosi2_out_p2_s" + str(s) + "_N" + str(N) + "_" + dom.strip() + "_rep" + str(rep) + "_" + args.suf +  ".txt"

				#file that will store output of ms_summarize.R
				summary_outfile = args.o + "cosi2_output/cosi2_summary_p2_s" + str(s) + "_N" + str(N) + "_n" + str(samp_sizes[0]) + "_" + dom.strip() + "_L" + args.L + "_rep" + str(rep) + "_" + args.suf + ".txt"

				R_args = [args.X, outfile, summary_outfile.replace(".txt",".tmp.txt"), args.L, args.w, args.W, str(num_pops), selPop]
				for n in samp_sizes: 
					R_args.append(str(int(n) * 2))

				awk_args = ['awk \'{printf \"2', str(rep), str(s), str(N), dom, '%s', ' '.join(samp_sizes), '\\n\",$0}\'', summary_outfile.replace(".txt",".tmp.txt"),'>',summary_outfile]
				cmd = "env COSI_NEWSIM=1 COSI_LOAD_TRAJ=" + traj_file + " coalescent -p " + param_file + " --nsims " + args.reps + " -m > " + outfile + " && " + " ".join(R_args) + " && " + " ".join(awk_args) + " && rm " + summary_outfile.replace(".txt",".tmp.txt")				
				print(cmd)
				#pp = subprocess.Popen(cmd,shell = True)



for s in tetTraj.s.unique():
	for N in tetTraj.N.unique():
		for dom in tetTraj.dom.unique():
			df=tetTraj.loc[(tetTraj.s == s) & (tetTraj.N == N) & (tetTraj.dom == dom),]
			for rep in df.rep.unique():

				# find generation at which selection stopped from the sweep pop_event and the population in which the sweep will occur
				events = args.events.split(",")
				sweeps = [s for s in events if "sweep" in s]
				assert len(sweeps) == 1, "Multiple sweeps specified in events.  Not currently supported."
				selStop = sweeps[0].split()[3]
				selPop = sweeps[0].split()[2]

				#Retrieve and format data for each rep of a set of run parameters
				# print("Running: ploidy = 2, s = " + str(s) + ", N = " + str(N) + ", Dominance = " + str(dom.strip()) + ", Rep = " + str(rep))
				df1=df.loc[df.rep == rep,("gen","freq")]
				selGen = len(df1.index) # Number of generations over which selection will occur
				df1['gen'] = abs(df1['gen'].astype(int) - selGen) + int(selStop)
				df1=df1.sort_values(by=["gen"])
				df1=df1.loc[(df1.freq.diff() < 0.15),] #any row that differs from previous row by more than 0.15 is likely an error...noticed several odd data points when visualizing data in R.
				if len(df1.loc[(df1.freq.diff() > 0.15),].index) > 3: # Should not find more than a couple artifactual data points
					df1['dif']= df1.freq.diff()
					print(df1)

				maxFreq=df1['freq'].max()

				# Write parameter file
				param_file = args.o + "cosi2_input/param_p4_s" + str(s) + "_N" + str(N) + "_" + dom.strip() + "_rep" + str(rep) + "_" + args.suf + ".txt"
				param_File = open(param_file,'w')
				param_File.write("# p = 4; s = " + str(s) + "; N = " + str(N) + "\n")
				param_File.write("length " + args.L + "\n\n")
				param_File.write("mutation_rate " + args.mu + "\n\n")
				param_File.write("recomb_file " + args.r + "\n\n")

				#Define populations
				for i, samp in enumerate(samp_sizes):
					if str(i) == selPop:
						param_File.write("pop_define " + str(i + 1) + " sweep" + str(i + 1) + "\n")
					else:
						param_File.write("pop_define " + str(i + 1) + " nosweep" + str(i + 1) + "\n")
					if args.N != "-9":
						pop_sizes = args.N.split(",")
						assert len(pop_sizes) >= selPop + 1
						param_File.write("pop_size " + str(i + 1) + " " + str(int(pop_sizes[i]) * 4) + "\n")
					
					##### IS THIS VALID? MULTIPLYING POPULATION SIZE BY PLOIDY? ####
					else: 
						param_File.write("pop_size " + str(i + 1) + " " + str(N * 4) + "\n")
					##### NEEDS TO BE CHECKED #####

					param_File.write("sample_size " + str(i + 1) + " " + str(int(samp * 4)) + "\n\n")

				#Define population events
				for event in events:
					event = event.split()
					if event[0] == "sweep":
						event[3] = str(selGen + int(selStop))
						event.append(str(maxFreq))
					elif event[0] == "split":
						T = int(event[4][1:])
						if event[4].startswith("B"):
							event[4] = str(T + selGen + int(selStop))
						elif event[4].startswith("D"):	
							event[4] = str(int(selStop) + T)
						elif event[4].startswith("A"):
							assert int(selStop) > T, "Negative (future) time for event: %r" % event[1]
							event[4] = str(int(selStop) - T)
					else:
						print("Event of type '" + event[0] + "' is not currently implemented.")
					param_File.write("pop_event " + " ".join(event) + "\n\n")
				param_File.close()

				#This files will hold the sweep trajectory for input into cosi2
				traj_file = args.o + "cosi2_input/traj_p4_s" + str(s) + "_N" + str(N) + "_" + dom.strip() + "_rep" + str(rep) + "_" + args.suf + ".txt"
				df1.to_csv(traj_file, index=False,header=None, sep=" ")

				#File that will store output of cosi2
				outfile = args.o + "cosi2_output/cosi2_out_p4_s" + str(s) + "_N" + str(N) + "_" + dom.strip() + "_rep" + str(rep) + ".txt"

				#file that will store output of ms_summarize.R
				summary_outfile = args.o + "cosi2_output/cosi2_summary_p4_s" + str(s) + "_N" + str(N) + "_n" + str(samp_sizes[0]) + "_" + dom.strip() + "_L" + args.L + "_rep" + str(rep) + "_" + args.suf + ".txt"

				# Prepare R command
				R_args = [args.X, outfile, summary_outfile.replace(".txt",".tmp.txt"), args.L, args.w, args.W, str(num_pops), selPop]
				for n in samp_sizes: 
					R_args.append(str(int(n) * 4))

				# Prepare awk command to add traj info to R output
				awk_args = ['awk \'{printf \"4', str(rep), str(s), str(N), dom, '%s', ' '.join(samp_sizes), '\\n\",$0}\'', summary_outfile.replace(".txt",".tmp.txt"), '>',summary_outfile]

				# Build master command
				cmd = "env COSI_NEWSIM=1 COSI_LOAD_TRAJ=" + traj_file + " coalescent -p " + param_file + " --nsims " + args.reps + " -m > " + outfile + " && " + " ".join(R_args) + " && " + " ".join(awk_args) + " && rm " + summary_outfile.replace(".txt",".tmp.txt")

				print(cmd)
				#pp = subprocess.Popen(cmd,shell = True)
