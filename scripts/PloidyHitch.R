#title: "PloidyHitch.R"
#author: "Patrick Monnahan"
#date: "11/23/2018"


# Load required packages
library(assertthat)
library(dplyr)
library(tidyr)
library(PopGenome)
library(stringr)
library(data.table)
library(magrittr)
library(doParallel)
library(parallel)

# Read in arguments
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

outdir = args[1]
outfile = args[2]
path_to_mssel = args[3]
dom_idx = args[4]
num_tries = args[5]
current_df = args[6] # file containing data of most current df of results.  used to keep track of how many reps we've done of each param set
num_cores = detectCores()

cl <- makeCluster(num_cores - 1) # Do we need to leave one core to handle the results that are being returned??
registerDoParallel(cl)

##### BEGIN: Set parameter values to explore #####
Ploidy = c(2, 4, 8)
pop_size = c(1000, 10000) # Keep this value above 100 and below 1000000 (computation time will increase with increasing pop_size)
selection_coeff = c(0.1, 0.01, 0.001) # Keep this between 0 and 1

dominance = c(0.1, 0.4, 0.5, 0.6, 0.9, 1, 0) #
dominance = dominance[dom_idx]

seq_len = 1000000 # Length of sequence that we will simulate with mssel.  Increasing this value will increase computation time.
mutation_rate = c(1e-8, 1e-7) # per-base mutation rate; mu
recomb_rate = c(1e-8, 1e-7) # per-base recombination rate; r
sampGen = c(1, 1000, 10000)
fuseTime = c(1, 1000, 10000)
samp_num = 10 # number of individuals to sample
num_reps = 5
max_gens = 9999
end_freq = 0.99

##### END: Set parameter values to explore #####

####### Begin: Define Functions ##########

# Retrieves dominance coefficients given ploidy and a single dominance scalar
getDomCoefs = function(dom_scalar, ploidy = 100){
  if (dom_scalar == 0){
    coeffs = rep(0, ploidy - 1)
  }
  else if (dom_scalar == 1){
    coeffs = rep(1, ploidy - 1)
  }
  else {
    dom = (1 - dom_scalar) - 0.5 
    dom = abs((dom * 10 + sign(dom)) ^ sign(dom))
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
      cat(p)
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
  mig_str = paste(npop, "0", n * p, n * p, "0 -ej", fuse_time, "2 1")
  
  # Prepare input/output
  writeTraj(traj_file, trajectory, npop, selPop, time_scale, name, sampleGen)
  
  # Format argument strings
  args1 = paste((2 * p * n), 1, n * p, n * p, traj_file, L * selPos, "-r", rho, L, "-t", theta, "-I", mig_str, ">", outfile)
  
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
      df$snp.start = as.numeric(str_split_fixed(rownames(df)," ",4)[,2])
      df$snp.end = as.numeric(str_split_fixed(str_split_fixed(rownames(df)," ",4)[,4]," ",2)[,1])
      df$rep1 = str_split_fixed(str_split_fixed(rownames(df),"_",3)[,3], "[.]",2)[,1]
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
  sim = readMS(in_file)
  positions = read.ms.output2(file.ms.output = in_file)
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
  sim.slide = sliding.window.transform(sim, width = winSize, jump = winSize * slideRate, type = 1, whole.data = FALSE)
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
  print("Done")
  
  return(df)
}

####### End: Define Useful Functions ##########

####### Begin: simulate ########

# Generate data table containing parameter sets of interest
params = expand.grid(traj_ploidy = Ploidy, sim_ploidy = Ploidy, s = selection_coeff, recomb = recomb_rate, N = pop_size, mu = mutation_rate, sampGen = sampGen, dom = dominance, fuseGen = fuseTime)
params %<>% dplyr::slice(rep(row_number(), num_reps))

params %<>% group_by(traj_ploidy, sim_ploidy, s, recomb, sampGen, N, mu, sampGen, dom, fuseGen) %>% mutate(rep = 1:n())

if (cdf != -9){
cdf = read.table(current_df, head = T)

done = cdf %>% select(s,dom, recomb, N, mu, sim_ploidy , traj_ploidy, sampGen, fuseGen, rep) %>% distinct() 

params %<>% anti_join(., done)

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

  # Get sweep trajectory
  new_traj = getTraj(s, dom, traj_ploidy, 0.05, N, end_freq, max_gens, num_tries)
  
  ms_outfile = paste(outdir, "/tp", traj_ploidy, "_sp", sim_ploidy, "_do", dom, "_se", s, "_NN", N, "_mu", mu, "_re", r, "_sG", sampGen, "_fG", fuseGen, "_rep", rep,  "_msel.out", sep = "")
  fin_outfile = paste(outdir, "/tp", traj_ploidy, "_sp", sim_ploidy, "_do", dom, "_se", s, "_NN", N, "_mu", mu, "_re", r, "_sG", sampGen, "_fG", fuseGen, "_rep", rep, "_smry.txt", sep = "")
  
  if (new_traj == -9){
    # print(paste("Beneficial allele was lost due to drift for", num_tries, "consecutive attempts"))
    print(params[j,])
  } else if (new_traj == -8){
    # print(paste("Beneficial allele did not fix before the maximum number of generations (", max_gens, ")."))
    print(params[j,])
  } else {
    new_traj = new_traj$freq
    
    infile = msselRun(N = N, n = samp_num, trajectory = new_traj, outfile = ms_outfile, L = seq_len, mu = , r = , ploidy = sim_ploidy, ms = path_to_mssel, sampleGen = sampGen, fuseGen = fuseGen)
  
    # calculate population genetic metrics in sliding windows across simulated region
    ndat = msselCalc(infile, numWindows = N / 50, rep(sim_ploidy * samp_num, 2), Nsites = seq_len)
    ndat %<>% mutate(s = s, dom = dom, recomb = r, N = N, mu = mu, sim_ploidy = sim_ploidy, traj_ploidy = traj_ploidy, sampGen = sampGen, fuseGen = fuseGen, rep = rep)
    
    # write.table(ndat, fin_outfile, row.names = F, quote = F)
  }
  return(ndat)
}

# Call to parallel
dat = mclapply(1:nrow(params), FUN = mssel_parallel, mc.cores= num_cores)

write.table(dat, paste(outdir, outfile), row.names=F, quote=F)
####### END: simulate ########



