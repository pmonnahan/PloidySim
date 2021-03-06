#title: "PloidyHitch.R"
#author: "Patrick Monnahan"
#date: "11/23/2018"


# Load required packages
library(assertthat)
library(plyr)
library(dplyr)
library(tidyr)
library(PopGenome)
library(stringr)
library(data.table)
library(magrittr)
library(doParallel)
library(parallel)
library(inline)
library(wrapr)

# Read in arguments
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

outdir = args[1]
outfile = args[2]
path_to_mssel = args[3]
dom_idx = as.numeric(args[4])
num_tries = as.numeric(args[5]) # Number of attempts to simulate the allele frequency trajectory.  If end_freq is not reached within max_gens, then another 'try' is made.
current_df = args[6] # file containing data of most current df of results.  used to keep track of how many reps we've done of each param set
recomb_rate = as.numeric(args[7]) # per-base recombination rate; r
num_reps = as.numeric(args[8])
num_drift_tries = 10000
num_cores = detectCores() - 1

#make tmpdir
tmpdir = paste(outdir,"/tmp", round(runif(1, 5000,1000000000)), sep="")
dir.create(tmpdir, showWarnings = FALSE)
setwd(tmpdir)

# cl <- makeCluster(num_cores) # Do we need to leave one core to handle the results that are being returned??
# registerDoParallel(cl)

##### BEGIN: Set parameter values to explore #####
Ploidy = c(2, 4, 8)
pop_size = c(1000, 10000) # Keep this value above 100 and below 1000000 (computation time will increase with increasing pop_size)
# selection_coeff = c(0.1, 0.01) # Keep this between 0 and 1
selection_coeff = c(0.1, 0.01) # Keep this between 0 and 1

# drift_gen = c(0, 1000, 10000)
drift_gen = c(0)

dominance = c(0.1, 0.5, 0.9, 1, 0) #
dominance = dominance[dom_idx]

seq_len = 1000000 # Length of sequence that we will simulate with mssel.  Increasing this value will increase computation time.
mutation_rate = c(1e-8) # per-base mutation rate; mu
sampGen = c(1, 1000, 10000)
fuseTime = c(1, 1000, 10000)
samp_num = as.numeric(args[9]) # number of individuals to sample
samp_alleles = 20
max_gens = 99999 # Number of generations to allow for simulation of allele frequency trajectory.  If end_freq is not reached before this point, then simulation is aborted.
end_freq = 0.99

same_ploidy = TRUE #Set to TRUE if you want traj_ploidy == sim_ploidy

##### END: Set parameter values to explore #####

####### Begin: Define Functions ##########
read.ms.output <- function( txt=NA, file.ms.output=NA,MSMS=FALSE) {
  
  if( !is.na(file.ms.output) ) txt <- scan(file=file.ms.output,
                                           what=character(0), sep="\n", quiet=TRUE)
  if( is.na(txt[1]) ){
    print("Usage: read.ms.output(txt), or read.ms.output(file=filename)")
    return()
  }
  
  
  if(MSMS[1]==FALSE){
    nsam   <- as.integer(strsplit(txt[1], split=" ")[[1]][2] )
    ndraws <- as.integer(strsplit(txt[1], split=" ")[[1]][3] )
  }
  
  #print(strsplit(txt[1], split=" "))
  
  h         <- numeric()
  result    <- list()
  gamlist   <- list()
  positions <- list()
  
  #marker <- grep("prob",txt)
  #probs <- sapply(strsplit(txt[marker], split=":"), function(vec) as.numeric(vec[2]))
  #marker <- grep("time",txt)
  #times <- sapply(strsplit(txt[marker], split="\t"), function(vec){ as.numeric(vec[2:3])} )
  times <- NaN
  probs <- NaN
  
  ## THE OUTPUT TEXT FOR EACH DRAW SHOULD CONTAIN THE WORD "segsites"
  marker <- grep("segsites", txt)
  
  if(MSMS[1]!=FALSE){ndraws <- length(marker);nsam <- MSMS$nsam} # MSMS
  stopifnot(length(marker) == ndraws)
  ## GET NUMBERS OF SEGREGATING SITES IN EACH DRAW
  segsites <- sapply(strsplit(txt[marker], split=" "), function(vec) as.integer(vec[2]) )
  
  for(draw in seq(along=marker)) {
    # if(!(draw %% 100)) cat(draw, " ")
    if(segsites[draw] > 0) {
      tpos <- strsplit(txt[marker[draw]+1], split=" ")
      positions[[draw]] <- as.numeric( tpos[[1]][ 2:(segsites[draw]+1) ] )
      
      haplotypes <- txt[(marker[draw] + 2):(marker[draw] + 2 + nsam - 1)]
      
      haplotypes <- strsplit(haplotypes, split="")
      
      h <- sapply(haplotypes, function(el) c(as.integer(el)))
      
      ## IF THERE'S 1 SEGREGATING SITE, THIS WON'T BE A MATRIX 
      
      if(segsites[draw] == 1) h <- as.matrix(h)
      ## OTHERWISE, IT NEEDS TO BE TRANSPOSED
      else h <- t(h)
    }
    else {
      h <- matrix(nrow=nsam, ncol=0)
      positions[[draw]] <- NA	
    }
    
    gamlist[[draw]] <- h
    stopifnot(all(dim(h) == c(nsam, segsites[draw]))) 
  }
  
  list(segsites=segsites, gametes=gamlist, probs=probs, times=t(times), positions=positions, nsam=nsam, nreps=ndraws ) 
}

readMS2 <- function(file, big.data=FALSE){
  dirname =  paste("SwapMS", round(runif(1, 5000,1000000000)), sep="")
  dir.create(dirname)
  if(!big.data){
    out     <- read.ms.output(file.ms.output=file)
    gametes <- out$gametes
    for(xx in 1:length(gametes)){
      
      d <- gametes[[xx]]
      d <- list(matrix=d,positions=NaN)
      samplename <- paste("ms_sample_",xx,".RD",sep="")
      save(d,file= file.path (dirname,samplename) )
      
    }
    test <- readData(dirname, SNP.DATA=F, FAST=TRUE, format="RData", big.data=big.data)
    unlink(dirname,recursive=F)
    return(test)
  }# end of !big.data
  
  if(big.data){
    out     <- read.big.ms.output(file)
    gametes <- out$gametes
    
    for(xx in 1:length(gametes)){
      open(gametes[[xx]])
      d <- gametes[[xx]][,]
      close(gametes[[xx]])
      d <- list(matrix=d,positions=NaN)
      samplename <- paste("ms_sample_",xx,".RD",sep="")
      save(d,file= file.path (dirname,samplename) )
      
    }
    test <- readData(dirname, SNP.DATA=F, FAST=TRUE, format="RData", big.data=big.data)
    unlink(dirname,recursive=F)
    return(test)
  }}#end of big.data
  
  
withErrorTracing = function(expr, silentSuccess=FALSE) {
  hasFailed = FALSE
  messages = list()
  warnings = list()
  
  errorTracer = function(obj) {
    
    # Storing the call stack 
    calls = sys.calls()
    calls = calls[1:length(calls)-1]
    # Keeping the calls only
    trace = limitedLabels(c(calls, attr(obj, "calls")))
    
    # Printing the 2nd and 3rd traces that contain the line where the error occured
    # This is the part you might want to edit to suit your needs
    print(paste0("Error occuring: ", trace[length(trace):1][2:3]))
    
    # Muffle any redundant output of the same message
    optionalRestart = function(r) { res = findRestart(r); if (!is.null(res)) invokeRestart(res) }
    optionalRestart("muffleMessage")
    optionalRestart("muffleWarning")
  }
  
  vexpr = withCallingHandlers(withVisible(expr),  error=errorTracer)
  if (silentSuccess && !hasFailed) {
    cat(paste(warnings, collapse=""))
  }
  if (vexpr$visible) vexpr$value else invisible(vexpr$value)
}

# Retrieves dominance coefficients given ploidy and a single dominance scalar
getDomCoefs = function(dom_scalar, ploidy = 100){
  if (dom_scalar == 0){
    coeffs = rep(0, ploidy - 1)
  }
  else if (dom_scalar == 1){
    coeffs = rep(1, ploidy - 1)
  }
  else {
    dom = 0.5 - dom_scalar # Center around additivity (0.5)
    dom = abs((dom * 10) ^ sign(dom))
    coeffs = (seq(1, ploidy - 1) / ploidy) ^ dom
  }
  return(coeffs)
}

# get trajectory and associated information given ploidy, dom_scalar, and selection coefficient
getTraj = function(s, dom_scalar, ploidy, start_freq, N = -9, end_freq = 0.99, maxGens = 9999, maxTries = 10){
  dom_coeffs = getDomCoefs(dom_scalar, ploidy)
  attempts = 0
  p = start_freq
  k = ploidy
  while (attempts <= maxTries & p < end_freq){
    df = data.frame("s" = s, "dom" = dom_scalar, "gen" = 0, "freq" = start_freq, "dp1" = 0, "w.bar" = 0, "ploidy" = ploidy)
    p = start_freq
    dp1 = 0
    gen = 0
    het_fits = 1 + dom_coeffs * s
    fits = c(1, het_fits, 1 + s) / (1 + s)
    while (p < end_freq & p > 0 & gen <= maxGens){
      q = 1 - p
      i = seq(0, k)
      Num = sum(choose(k, i) * ((k - i) / k) * (p ^ (k - i)) * (q ^ i) * rev(fits))
      Den = sum(choose(k, i) * (p ^ (k - i)) * (q ^ i) * rev(fits))
      p_prime = Num / Den
      if (N != -9){
        p_prime = sum(rbinom(N, k, p_prime)) / (k * N) 
      }
      df = rbind(df, c(s, dom_scalar, gen, p, dp1, Den, k))
      dp1 = p_prime - p
      p = p_prime
      gen = gen + 1
    }
    df = rbind(df, c(s, dom_scalar, gen, p, dp1, Den, k))
    attempts = attempts + 1
  }
  if (attempts > maxTries & nrow(df) <= maxGens){
    print(paste("Beneficial allele was lost due to drift for", maxTries, "consecutive attempts"))
    return(-9)
  } else if (nrow(df) > maxGens){
    print(paste("Beneficial allele did not fix before the maximum number of generations (", maxGens, ")."))
    return(-8)
  } else{
    return(df[-1,])
  }
}

Drift = function(N, ploidy, end, tries){
  drift = function(N, ploidy, end){
    C = N * ploidy # # of chromosomes
    p = 0 # Current frequency
    traj = c(1/C)
    while(p == 0){
      p = rbinom(1, C, (1 / C)) / (C) # Draw # of mutant alleles for next generation and divive by C to get freq
    }
    t = 1
    traj = c(traj, p)
    while(p != 0 & p != 1 & t < end){ # Loop till absorption or beginning of selection
      p = rbinom(1, C, p) / C
      t = t + 1
      traj = c(traj, p)
    }
    return(traj)
  }
  
  tt = 1
  if(end > 0){
    traj = drift(N, ploidy, end)
    # if (length(traj) == 0){traj = c(0)}
    # What is happening after tt exceeds tries????  Is this getting auto set to 0? such that start_freq is 0 for selection??
    while(tt < tries & (traj[length(traj)] == 0 | traj[length(traj)] == 1)){
      traj = drift(N, ploidy, end)
      tt = tt + 1
    }
  } else{
    traj = 1 / (N * ploidy)
  }
  return(traj)
}

# Format and write allele frequency trajectory to file for input into mssel
writeTraj = function(file, traj, numPops, selPops, timeScale, trajName = "rep1", startGen = 1){
  assert_that(length(traj) > 10, msg = "Error with trajectory")
  gen = seq(startGen, startGen + length(traj)) / timeScale
  fileConn = file(file)
  writeLines(c(paste("ntraj:", "1"), paste("npop:", numPops), paste("n:", length(traj) + 1, trajName)), con = fileConn)
  dat = data.frame("gen" = gen, "traj" = c(rev(traj), 0), "anc" = rep(0, length(traj) + 1))
  write.table(dat, file, quote = F, col.names = F, row.names = F, sep = "\t", append = T)
  close(fileConn)
}

# Run mssel
msselRun = function(N, n, trajectory, outfile = -9, L = 1000000, mu = 1e-8, r = 1e-8, ploidy = 2, selPos = 0.5, npop = 2, selPop = 1, fuseGen = 1, sampleGen = 1, ms = "/Users/pmonnahan/Documents/Research/code/dmc/mssel_modified/mssel"){ 
  options(scipen=999)
  if (outfile == -9){
    name = paste(getwd(), "/mssel.", sample(1:99999999, 1), sep = "")
    outfile = paste(name, ".ms.out", sep = "")
    traj_file = paste(name, ".traj.txt", sep = "")
  } else {
    name = tools::file_path_sans_ext(outfile)
    traj_file = paste(name, ".traj.txt", sep = "")
  }
  p = ploidy
  time_scale = 2 * N * p
  theta = 2 * N * p * mu * L 
  rho = 2 * N * p * r * L
  fuse_time =  (length(trajectory) + sampleGen + fuseGen) / time_scale
  mig_str = paste(npop, "0", n, n, "0 -ej", fuse_time, "2 1")
  
  # Prepare input/output
  writeTraj(traj_file, trajectory, npop, selPop, time_scale, name, sampleGen)
  
  # Format argument strings
  args1 = paste((2 * n), 1, n, n, traj_file, L * selPos, "-r", rho, L, "-t", theta, "-I", mig_str, ">", outfile)
  
  # Run mssel
  print(paste(ms, args1))
  cmd1 = system2(ms, args1)
  return(outfile)
}

# Use PopGenome package to calculate a number of popgen stats
msselCalc <- function(in_file, numWindows, samp_sizes, Nsites, outgroup = 21, slideRate=0.5, selPop=1, linkage_stats = c("Kelly.Z_nS"), neutrality_stats = c("Tajima.D", "Fay.Wu.H", "Zeng.E")){
  # Define some necessary functions
  read.ms.output2 <- function(txt=NA, file.ms.output=NA, MSMS=FALSE) {
    
    if( !is.na(file.ms.output) ) txt <- scan(file=file.ms.output,
                                             what=character(0), sep="\n", quiet=TRUE)
    if( is.na(txt[1]) ){
      print("Usage: read.ms.output(txt), or read.ms.output(file=filename)")
      return()
    }
    
    
    if(MSMS[1]==FALSE){
      nsam   <- as.integer(strsplit(txt[1], split=" ")[[1]][2] )
      ndraws <- as.integer(strsplit(txt[1], split=" ")[[1]][3] )
    }
    
    #print(strsplit(txt[1], split=" "))
    
    h         <- numeric()
    result    <- list()
    gamlist   <- list()
    positions <- list()
    
    #marker <- grep("prob",txt)
    #probs <- sapply(strsplit(txt[marker], split=":"), function(vec) as.numeric(vec[2]))
    #marker <- grep("time",txt)
    #times <- sapply(strsplit(txt[marker], split="\t"), function(vec){ as.numeric(vec[2:3])} )
    times <- NaN
    probs <- NaN
    
    ## THE OUTPUT TEXT FOR EACH DRAW SHOULD CONTAIN THE WORD "segsites"
    marker <- grep("segsites", txt)
    
    if(MSMS[1]!=FALSE){ndraws <- length(marker);nsam <- MSMS$nsam} # MSMS
    
    
    stopifnot(length(marker) == ndraws)
    
    ## GET NUMBERS OF SEGREGATING SITES IN EACH DRAW
    segsites <- sapply(strsplit(txt[marker], split=" "), function(vec) as.integer(vec[2]) )
    
    for(draw in seq(along=marker)) {
      # if(!(draw %% 100)) cat(draw, " ")
      if(segsites[draw] > 0) {
        tpos <- strsplit(txt[marker[draw]+1], split=" ")
        positions[[draw]] <- as.numeric( tpos[[1]][ 2:(segsites[draw]+1) ] )
        
        haplotypes <- txt[(marker[draw] + 2):(marker[draw] + 2 + nsam - 1)]
        
        haplotypes <- strsplit(haplotypes, split="")
        
        h <- sapply(haplotypes, function(el) c(as.integer(el)))
        
        ## IF THERE'S 1 SEGREGATING SITE, THIS WON'T BE A MATRIX 
        
        if(segsites[draw] == 1) h <- as.matrix(h)
        ## OTHERWISE, IT NEEDS TO BE TRANSPOSED
        else h <- t(h)
        
        
      }
      else {
        h <- matrix(nrow=nsam, ncol=0)
        positions[[draw]] <- NA	
      }
      
      gamlist[[draw]] <- h
      stopifnot(all(dim(h) == c(nsam, segsites[draw]))) 
    }
    
    list(segsites=segsites, gametes=gamlist, probs=probs, times=t(times), positions=positions, nsam=nsam, nreps=ndraws ) 
  }
  makePops <- function(sampSize_list){
    popList = list()
    start = 1
    running_tot = 0
    for (i in 1:length(sampSize_list)){
      running_tot = running_tot + sampSize_list[i]
      popList[[i]] = seq(start,running_tot)
      start = start + sampSize_list[i]
    }
    return(popList)
  }
  munge <- function(inp, Nsites, positions, selPop, numWindows){
    coords = c("rep1","snp.start", "snp.end")
    DF = data.frame()
    parseInfo = function(info){
      df = as.data.frame(info)
      df$snp.start = as.numeric(str_split_fixed(rownames(df)," ",4)[,1])
      df$snp.end = as.numeric(str_split_fixed(rownames(df)," ",4)[,3])
      df$rep1 = 1
      return(df)
    }
    for (i in 1:length(inp@populations)){
      div = parseInfo(get.diversity(inp)[[i]])
      ld = parseInfo(get.linkage(inp)[[i]])
      neutrality_tests = parseInfo(get.neutrality(inp)[[i]])
      df = merge(div, ld, by = coords)
      df = merge(df, neutrality_tests[, c(coords, neutrality_stats)], by = coords)
      df$pop = as.character(i)
      DF=rbind(DF,df)
    }
    
    vars = c(coords, "nuc.diversity.within", linkage_stats, neutrality_stats, "pop")
    DF = melt(DF[,vars], id.vars = c(coords, "pop"))
    DF = dcast(DF, rep1 + snp.start + snp.end ~ variable + pop, value.var = "value")
    pos = melt(positions)
    colnames(pos)[2] = "rep1"
    pos$rep1=as.character(pos$rep1)
    ends = pos %>% group_by(rep1) %>% mutate(snp.end=row_number()) %>% dplyr::slice(unique(DF$snp.end)) %>% as.data.frame()
    starts = pos %>% group_by(rep1) %>% mutate(snp.start=row_number()) %>% dplyr::slice(unique(DF$snp.start)) %>% as.data.frame()
    colnames(ends)[1] = "bp.end"
    colnames(starts)[1] = "bp.start"
    ends$bp.end = ends$bp.end * Nsites
    starts$bp.start = starts$bp.start * Nsites
    DF = merge(DF,starts,by=c("snp.start","rep1"))
    DF = merge(DF,ends,by=c("snp.end","rep1"))
    
    DF = DF[order(DF$rep1, DF$snp.start),]
    
    dxy = dcast(melt(inp@nuc.diversity.between), Var2 ~ Var1, value.var="value")
    colnames(dxy)[2] = 'dxy'
    fst = dcast(melt(inp@nuc.F_ST.pairwise), Var2 ~ Var1, value.var="value")
    colnames(fst)[2] = 'fst'
    DF = cbind(DF, dxy)[,-2]
    DF = cbind(DF, fst)[,-2]
    DF$selPop = selPop
    return(DF)
  }
  
  # Read output of mssel simulation
  print("Reading mssel input...")
  sim = readMS2(in_file)
  print("Got polymorphism data...")
  positions = read.ms.output2(file.ms.output = in_file)
  print("Got position data...")
  positions = positions$positions
  winSize = round(length(positions[[1]]) / numWindows, digits = -1)
  print("Done")
  
  # Define populations...number of anc and der specified in command line?
  pops = makePops(samp_sizes)
  sim=set.populations(sim,pops) #Only works for two populations currently
  sim=diversity.stats(sim)
  sim=diversity.stats.between(sim)
  
  ## EXPERIMENTAL ##
  sim = set.outgroup(sim, new.outgroup = outgroup)
  ## ##########
  
  # Divide simulated data into windows (type=1: based on SNP counts; type=2: based on nucleotide counts)
  cat(paste("","Creating windows...", "", sep = "\n"))
  sim.slide = sliding.window.transform(sim, width = winSize, jump = winSize * slideRate, type = 1, whole.data = T) #whole.data must be set to F if using the ff package (i.e. if SNP.DATA is set to T)
  cat(paste("Done", "Calculating metrics...", "", sep = "\n"))
  sim.slide = diversity.stats(sim.slide)
  sim.slide = diversity.stats.between(sim.slide)
  sim.slide = linkage.stats(sim.slide)
  sim.slide = F_ST.stats(sim.slide)
  sim.slide = neutrality.stats(sim.slide, detail=TRUE)
  cat(paste("Done", "Munging data...", "", sep = "\n"))
  df = munge(sim.slide, Nsites, positions, selPop)
  df$mid = (df$bp.end + df$bp.start)/2
  df$bp.len = (df$bp.end - df$bp.start)
  df$Pi.1 = df$nuc.diversity.within_1 / df$bp.len
  df$Pi.2 = df$nuc.diversity.within_2 / df$bp.len
  df$dxy.1.2 = df[['dxy']] / df$bp.len
  df %<>% select(-one_of(c("Var2, Var2.1"))) %>% as.data.frame()
  print("Done")
  
  return(df)
}


####### End: Define Useful Functions ##########

####### Begin: simulate ########

# Generate data table containing parameter sets of interest
params = expand.grid(traj_ploidy = Ploidy, sim_ploidy = Ploidy, s = selection_coeff, recomb = recomb_rate, N = pop_size, mu = mutation_rate, sampGen = sampGen, dom = dominance, fuseGen = fuseTime, driftGen = drift_gen)
params %<>% dplyr::slice(rep(row_number(), num_reps))

params %<>% group_by(traj_ploidy, sim_ploidy, s, recomb, sampGen, N, mu, sampGen, dom, fuseGen, driftGen) %>% mutate(rep = 1:n())

if (same_ploidy){
  params %<>% filter(traj_ploidy == sim_ploidy)
  }

if (current_df != -9){
cdf = read.table(current_df, head = T)

done = cdf %>% select(s,dom, recomb, N, mu, sim_ploidy , traj_ploidy, sampGen, fuseGen, rep, driftGen) %>% distinct() 

params %<>% anti_join(., done)

remove(cdf)
remove(done)

}

# Run mssel in parallel for all parameter sets in params

mssel_parallel = function(j) {

  # parse parameter info
  traj_ploidy = params[j,]$traj_ploidy
  sim_ploidy = params[j,]$sim_ploidy
  dom = params[j,]$dom
  s = params[j,]$s
  N = params[j,]$N
  mu = params[j,]$mu
  r = params[j,]$recomb
  sampGen = params[j,]$sampGen
  fuseGen = params[j,]$fuseGen
  rep = params[j,]$rep
  driftGen = params[j,]$driftGen
  
  print(paste("hey", N, traj_ploidy, driftGen, num_drift_tries))
  drift_traj = Drift(N, traj_ploidy, driftGen, num_drift_tries)
  start_freq = drift_traj[length(drift_traj)]
  
  #Create filenames
  ms_outfile = paste(outdir, "/tp", traj_ploidy, "_sp", sim_ploidy, "_do", dom, "_se", s, "_NN", N, "_mu", mu, "_re", r, "_sG", sampGen, "_fG", fuseGen, "_dG", driftGen ,"_rep", rep,  "_msel.out", sep = "")
  fin_outfile = paste(outdir, "/tp", traj_ploidy, "_sp", sim_ploidy, "_do", dom, "_se", s, "_NN", N, "_mu", mu, "_re", r, "_sG", sampGen, "_fG", fuseGen, "_dG", driftGen, "_rep", rep, "_smry.txt", sep = "")
  
  # Get sweep trajectory
  sel_traj = tryCatch({
    getTraj(s, dom, traj_ploidy, start_freq, N, end_freq, max_gens, num_tries)
  }, 
  warning = function(w){
    message("Warning: getTraj")
    message(params[j,])
    return(NULL)
    },
  error = function(e){
    message("Error: getTraj")
    message(params[j,])
    return(NULL)
    })
  
  # Run mssel
  if (is.data.frame(sel_traj)){
    infile = tryCatch({
    sel_traj = sel_traj$freq
    new_traj = c(drift_traj, sel_traj[2:length(sel_traj)])
    if (samp_num != -9){
      msselRun(N = N, n = samp_num * sim_ploidy, trajectory = new_traj, outfile = ms_outfile, L = seq_len, mu = mu, r = r, ploidy = sim_ploidy, ms = path_to_mssel, sampleGen = sampGen, fuseGen = fuseGen)
    }
    else{
      msselRun(N = N, n = samp_alleles, trajectory = new_traj, outfile = ms_outfile, L = seq_len, mu = mu, r = r, ploidy = sim_ploidy, ms = path_to_mssel, sampleGen = sampGen, fuseGen = fuseGen)
    }
    }, 
    warning = function(w){
      message("Warning: msselRun")
      message(params[j,])
      return(NULL)
    },
    error = function(e){
      message("Error: msselRun")
      message(params[j,])
      return(NULL)
    })} else {
      message("Trajectory failure")
      message(params[j,])
    }
  # calculate population genetic metrics in sliding windows across simulated region
  if (exists("infile") && !is.null(infile)){
    ndat = tryCatch({
      if (samp_num != -9){
    msselCalc(infile, numWindows = N / 50, rep(sim_ploidy * samp_num, 2), Nsites = seq_len)
      }
      else{
        msselCalc(infile, numWindows = N / 50, rep(samp_alleles, 2), Nsites = seq_len)
      }
      }, 
    error = function(e){
      message("Error: msselCalc")
      errorStorage <- capture.output(tryCatch({
        withErrorTracing({msselCalc(infile, numWindows = N / 50, rep(sim_ploidy * samp_num, 2), Nsites = seq_len)})
      }, error = function(e){
        e <<- e
        message("ERROR: ", e$message, "\nin ")
        message(e$call)
      }))
      return(NULL)
    })
  }
  
  if (!exists("ndat", inherits=F) || is.null(ndat)){
    ndat = params[j,]
  } else{
    ndat %<>% mutate(s = s, dom = dom, recomb = r, N = N, mu = mu, sim_ploidy = sim_ploidy, traj_ploidy = traj_ploidy, sampGen = sampGen, fuseGen = fuseGen, driftGen = driftGen, start_freq = start_freq, rep = rep) %>% select(-c(Var2, Var2.1))
    }
    return(ndat)
}

# Call to parallel

includes <- '#include <sys/wait.h>'
code <- 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
wait <- cfunction(body=code, includes=includes, convention='.C')

for (i in 1:ceiling((nrow(params)/num_cores))){
  
  if (i == ceiling((nrow(params)/num_cores))){
    end_idx = nrow(params)
    start_idx = (i - 1) * num_cores
  } else {
    end_idx = i * num_cores
    start_idx = (end_idx - num_cores) + 1
  }
  
  progress = round(i/(nrow(params)/num_cores), digits = 4) * 100
  print(paste0("Running jobs: ", start_idx, "-", end_idx, "; (", progress, "% complete)"))
  ndat = mclapply(start_idx:end_idx, FUN = mssel_parallel, mc.cores = num_cores)
  wait()
  # print(ndat)
  ndat = as.data.frame(do.call(rbind.fill, ndat))
  if (i==1){
    dat = ndat
    Dat = ndat
  } else if (i %% 100){
    write.table(dat, paste(outdir, outfile, "prm", i, ".txt", sep = ""), row.names=F, quote=F)
    dat = ndat
    Dat = rbind.fill(ndat, Dat)
  } else{
    dat = rbind.fill(ndat, dat)
    Dat = rbind.fill(ndat, Dat)
  }
}


write.table(dat, paste(outdir, outfile, "prm_fin", ".txt", sep = ""), row.names=F, quote=F)
write.table(Dat, paste(outdir, outfile, ".txt", sep = ""), row.names=F, quote=F)
unlink(tmpdir, recursive=T)
####### END: simulate ########



