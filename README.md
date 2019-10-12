# PloidySim

Coalescent simulations (using mssel) to compare population genomic footprints of selection across different ploidy levels.

Associated citation: \
The effect of autopolyploidy on population genetic signals of hard sweeps \
Patrick Monnahan, Yaniv Brandvain \
bioRxiv 753491; doi: https://doi.org/10.1101/753491

# Summary

The primary script used to generate allele frequency trajectories, run mssel (simulate polymorphism data), and calculate basic population genetic metrics (excluding haplotype metrics; accomplished with scripts/calcREHH.R) is scripts/PloidyHitch.R.

The mssel output files generated from PloidyHitch.R can be subsequently used for calcREHH.R, which will calculate the iHS and XPEHH haplotype-based statistics.

Figures in the above citation can be created with the script analysis/PaperFigures.Rmd.  The specific data files referenced in this script can be found at the following dryad repository:

Other useful code for visualizing the output of PloidyHitch.R and calcREHH.R can be found in analysis/mssel_analysis.Rmd.

There are a number of required R packages that are necessary to running the code contained in this repository.  They are:
* rehh
* PopGenome
* assertthat
* tidyverse
* stringr
* data.table
* magrittr
* doParallel
* parallel
* inline
* wrapr
* rio
* readr 

Many of these packages may have additional dependencies that the user would need to acquire as well.

## PloidyHitch.R
This wrapper loops over all factorial combinations of parameters (jointly specified in the command line and top of the script) and for each combination, it generates an allele frequency trajectory, runs mssel, and calculates basic population genetic metrics, using the R package PopGenome.

Runs are parallelized across parameter combinations, which may cause issues on some operating systems or cluster configurations.  The number of processors can be set to 1 to bypass parallelization if necessary. 

### Usage
    Rscript PloidyHitch.R <output_directory> <output_file> <path_to_mssel> <dom_idx> <num_tries> <done_df> <recomb_rate> <num_reps> <num_drift_tries> <samp_num>

## calcREHH.R
### Usage
