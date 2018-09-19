#!/bin/bash
#PBS -l mem=48gb,nodes=1:ppn=24,walltime=24:00:00
#PBS -A brandvai
#PBS -m abe
#PBS -M pmonnaha@umn.edu
#PBS -o /home/brandvai/pmonnaha/OandE/PloidySim.o
#PBS -e /home/brandvai/pmonnaha/OandE/PloidySim.e
#PBS -q mesabi

# Load modules
module load parallel
module load python3/3.4

echo -e "rep\tploidy\tstartfreq\tN\ts\tdom\tgen\tNgen\tfreq" > /home/brandvai/pmonnaha/PloidySim/Ploidy_Forward_Sim_results.txt

TASKS="/home/brandvai/pmonnaha/PloidySim/Ploidy_For_Sim_commands_cluster.txt"

cat ${PBS_NODEFILE}

# Uncomment the next line for running parallel across multiple nodes
#parallel -v -v -t --jobs 24 --sshloginfile ${PBS_NODEFILE} --workdir ${PWD} < ${TASKS}
# Uncomment the next line for single-node parallel
parallel < ${TASKS}