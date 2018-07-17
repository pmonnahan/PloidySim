#!/usr/bin/env Rscript

require("PopGenome")
require("stringr")
require("data.table")
require("dplyr")


args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
  # default window size
  args[3] = 50
}
  

# HARDCODED ARGS TO BE RELAXED LATER
Nanc=5 #Number of ancestral samples
Nder=5 #Number of derived/sweep samples
Nsites=50000 #Number of sites
winSize=50 #In terms of number of snps
slideRate = 0.5 # Proportion of winSize that we want to slide the window

# Necessary Functions
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


makePop.df=function(inp,Nsites,positions){
  DF = data.frame()
  
  for (i in 1:length(inp@populations)){
    df=as.data.frame(get.diversity(inp)[[i]])
    df$snp.start=as.numeric(str_split_fixed(rownames(df)," ",4)[,2])
    df$snp.end=as.numeric(str_split_fixed(str_split_fixed(rownames(df)," ",4)[,4]," ",2)[,1])
    df$rep=str_split_fixed(str_split_fixed(rownames(df),"_",3)[,3], "[.]",2)[,1]
    
    df1=as.data.frame(get.linkage(inp)[[i]])
    df1$snp.start=as.numeric(str_split_fixed(rownames(df1)," ",4)[,2])
    df1$snp.end=as.numeric(str_split_fixed(str_split_fixed(rownames(df1)," ",4)[,4]," ",2)[,1])
    df1$rep=str_split_fixed(str_split_fixed(rownames(df1),"_",3)[,3], "[.]",2)[,1]
    df = merge(df,df1,by=c("rep","snp.start","snp.end"))
    
    df$pop = as.character(i)
    DF=rbind(DF,df)
  }
  DF=dcast(setDT(DF), rep + snp.start + snp.end ~ pop,value.var = c("nuc.diversity.within","Kelly.Z_nS"))
  pos = melt(positions)
  colnames(pos)[2] = "rep"
  pos$rep=as.character(pos$rep)
  ends = pos %>% group_by(rep) %>% mutate(snp.end=row_number()) %>% slice(unique(DF$snp.end)) %>% as.data.frame()
  starts = pos %>% group_by(rep) %>% mutate(snp.start=row_number()) %>% slice(unique(DF$snp.start)) %>% as.data.frame()
  colnames(ends)[1] = "bp.end"
  colnames(starts)[1] = "bp.start"
  
  ends$bp.end=ends$bp.end * Nsites
  starts$bp.start=starts$bp.start * Nsites
  DF=merge(DF,starts,by=c("snp.start","rep"))
  DF=merge(DF,ends,by=c("snp.end","rep"))
  
  DF=DF[order(rep,snp.start)]
  dxy=dcast(melt(inp@nuc.diversity.between), Var2 ~ Var1, value.var="value")
  
  DF = cbind(DF,dxy)
  return(DF)
}

# Read output of cosi2 simulation
sim = readMS(args[1])
positions = read.ms.output(file.ms.output=args[1])
positions = positions$positions

# Define populations...number of anc and der specified in command line?

sim=set.populations(sim,list(c(as.character(seq(1,Nanc))),c(as.character(seq(Nanc+1,Nder + Nanc))))) #Only works for two populations currently
sim=diversity.stats(sim)
sim=diversity.stats.between(sim)

# Divide simulated data into windows
# type=1: Define windows based on SNP counts
# type=2: Define windows based on nucleotide counts
sim.slide = sliding.window.transform(sim,width=winSize, jump= winSize * slideRate,type=1,whole.data=FALSE)
sim.slide = diversity.stats(sim.slide)
sim.slide = diversity.stats.between(sim.slide)
sim.slide = linkage.stats(sim.slide)
sim.slide = F_ST.stats(sim.slide)

df = makePop.df(sim.slide, Nsites, positions)

write.table(df,args[2],row.names = FALSE, quote = FALSE)



