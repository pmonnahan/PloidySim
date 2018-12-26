#!/usr/bin/env Rscript

require("PopGenome")
require("stringr")
require("data.table")
require("dplyr")
require("ggplot2")

quickCosi1=function(cosi_output_file){
  co2=readMS(file=cosi_output_file)
  posb=read.ms.output(file.ms.output =cosi_output_file)
  posb2=posb$positions
  cow2=sliding.window.transform(co2,width=20,jump=10,type=1,whole.data = F)
  cow2 = diversity.stats(cow2)
  cow2 = diversity.stats.between(cow2)
  cow2 = linkage.stats(cow2)
  cow2 = F_ST.stats(cow2)
  df2 = makeOnePop.df(cow2,500000,posb2,1)
  ggplot(df2,aes(x=mid,y=Pi,color=rep)) + geom_line()
  return(df2)
}

quickCosiM=function(cosi_output_file){
  co2=readMS(file=cosi_output_file)
  posb=read.ms.output(file.ms.output =cosi_output_file)
  posb2=posb$positions
  cow2=sliding.window.transform(co2,width=20,jump=10,type=1,whole.data = F)
  cow2 = diversity.stats(cow2)
  cow2 = diversity.stats.between(cow2)
  cow2 = linkage.stats(cow2)
  cow2 = F_ST.stats(cow2)
  df2 = makeOnePop.df(cow2,500000,posb2,1)
  ggplot(df2,aes(x=mid,y=Pi,color=rep)) + geom_line()
  return(df2)
}

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

in_file = args[1]
out_file = args[2]
Nsites = as.numeric(args[3])
numWindows = as.numeric(args[4])
slideRate = as.numeric(args[5])
Npops = as.numeric(args[6])
selPop = as.numeric(args[7])
samp_sizes = as.numeric(c(args[8:length(args)]))

# in_file = "/Users/pmonnahan/Documents/Research/PloidySim/mssel_out_p2.0_s0.1_N1000000_additive_rep0Alt.txt"
# out_file = str_replace(in_file,"txt",".out")
# Nsites = 1000000
# numWindows = 20
# slideRate = 0.5
# Npops = 2
# selPop = 2
# samp_sizes = as.numeric(c(20,20))

# Necessary Functions
read.ms.output2 <- function( txt=NA, file.ms.output=NA,MSMS=FALSE) {
  
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

munge = function(inp,Nsites,positions,selPop, numWindows){
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
    df3 = as.data.frame(get.neutrality(inp)[[i]])
    df3$snp.start=as.numeric(str_split_fixed(rownames(df3)," ",4)[,2])
    df3$snp.end=as.numeric(str_split_fixed(str_split_fixed(rownames(df3)," ",4)[,4]," ",2)[,1])
    df3$rep=str_split_fixed(str_split_fixed(rownames(df3),"_",3)[,3], "[.]",2)[,1]
    
    df = merge(df, df3[, c("rep","snp.start", "snp.end", "Tajima.D")], by = c("rep", "snp.start", "snp.end"))
    df$pop = as.character(i)
    DF=rbind(DF,df)
  }
  
  vars = c("rep","snp.start","snp.end","nuc.diversity.within","Kelly.Z_nS","Tajima.D","pop")
  DF = melt(DF[,vars], id.vars = c("rep","snp.start", "snp.end", "pop"))
  DF = dcast(DF, rep + snp.start + snp.end ~ variable + pop, value.var = "value")
  pos = melt(positions)
  colnames(pos)[2] = "rep"
  pos$rep=as.character(pos$rep)
  ends = pos %>% group_by(rep) %>% mutate(snp.end=row_number()) %>% slice(unique(DF$snp.end)) %>% as.data.frame()
  starts = pos %>% group_by(rep) %>% mutate(snp.start=row_number()) %>% slice(unique(DF$snp.start)) %>% as.data.frame()
  colnames(ends)[1] = "bp.end"
  colnames(starts)[1] = "bp.start"
  ends$bp.end=ends$bp.end * Nsites
  starts$bp.start=starts$bp.start * Nsites
  DF = merge(DF,starts,by=c("snp.start","rep"))
  DF = merge(DF,ends,by=c("snp.end","rep"))
  
  DF = DF[order(DF$rep, DF$snp.start),]
  
  dxy = dcast(melt(inp@nuc.diversity.between), Var2 ~ Var1, value.var="value")
  colnames(dxy)[2] = 'dxy'
  fst = dcast(melt(inp@nuc.F_ST.pairwise), Var2 ~ Var1, value.var="value")
  colnames(fst)[2] = 'fst'
  print(str(dxy))
  print(str(fst))
  print(str(DF))
  DF = cbind(DF, dxy)
  DF = cbind(DF, fst)
  DF$selPop = selPop
  DF=DF[,c(2,3,1,10,11,4,5,6,7,8,9,11,13,15,16)]
  return(DF)
}

makeOnePop.df=function(inp,Nsites,positions,popIndex){
  df=as.data.frame(get.diversity(inp)[[popIndex]])
  df$snp.start=as.numeric(str_split_fixed(rownames(df)," ",4)[,2])
  df$snp.end=as.numeric(str_split_fixed(str_split_fixed(rownames(df)," ",4)[,4]," ",2)[,1])
  df$rep=str_split_fixed(str_split_fixed(rownames(df),"_",3)[,3], "[.]",2)[,1]
  df1=as.data.frame(get.linkage(inp)[[popIndex]])
  df1$snp.start=as.numeric(str_split_fixed(rownames(df1)," ",4)[,2])
  df1$snp.end=as.numeric(str_split_fixed(str_split_fixed(rownames(df1)," ",4)[,4]," ",2)[,1])
  df1$rep=str_split_fixed(str_split_fixed(rownames(df1),"_",3)[,3], "[.]",2)[,1]
  df = merge(df,df1,by=c("rep","snp.start","snp.end"))
  df$pop = as.character(popIndex)
  pos = melt(positions)
  print(head(pos))
  colnames(pos)[2] = "rep"
  pos$rep=as.character(pos$rep)
  print(head(unique(df$snp.end)))
  ends = pos %>% group_by(rep) %>% mutate(snp.end=row_number()) %>% slice(unique(df$snp.end)) %>% as.data.frame()
  starts = pos %>% group_by(rep) %>% mutate(snp.start=row_number()) %>% slice(unique(df$snp.start)) %>% as.data.frame()
  colnames(ends)[1] = "bp.end"
  colnames(starts)[1] = "bp.start"
  ends$bp.end=ends$bp.end * Nsites
  starts$bp.start=starts$bp.start * Nsites
  df=merge(df,starts,by=c("snp.start","rep"))
  df=merge(df,ends,by=c("snp.end","rep"))
  
  df=df[order(df$rep,df$snp.start),]
  df$rep = as.factor(df$rep)
  df$mid = (df$bp.end + df$bp.start)/2
  df$bp.len = (df$bp.end - df$bp.start)
  df$Pi = df$nuc.diversity.within / df$bp.len
  return(df)
}

# Read output of mssel simulation
sim = readMS(in_file)
positions = read.ms.output2(file.ms.output=in_file)
positions = positions$positions
winSize = round(length(positions[[1]]) / numWindows, digits = -1)

# Define populations...number of anc and der specified in command line?
pops = makePops(samp_sizes)
sim=set.populations(sim,pops) #Only works for two populations currently
sim=diversity.stats(sim)
sim=diversity.stats.between(sim)

## EXPERIMENTAL ##

sim = set.outgroup(sim, new.outgroup = seq(1,20))

## ##########

# Divide simulated data into windows
# type=1: Define windows based on SNP counts
# type=2: Define windows based on nucleotide counts

sim.slide = sliding.window.transform(sim, width = winSize, jump = winSize * slideRate, type = 1, whole.data = FALSE)
sim.slide = diversity.stats(sim.slide)
sim.slide = diversity.stats.between(sim.slide)
sim.slide = linkage.stats(sim.slide)
sim.slide = F_ST.stats(sim.slide)
# sim.slide = mult.linkage.stats(sim)
# sim.slide = neutrality.stats(sim.slide, new.outgroup = seq(0, samp_sizes[1]))
sim.slide = neutrality.stats(sim.slide, detail=TRUE)
# sim = detail.stats(sim)
# freq = sim@region.stats@minor.allele.freqs[[1]]
# freq.table <- list()
# freq.table[[1]] <- table(freq)
# sim.slisim <- splitting.data(sim.slide, positions= list(900:1100, 1200:2000), type=1)
# sim.slide = sweeps.stats(sim.slide)


df = munge(sim.slide, Nsites, positions, selPop)
df$rep = as.factor(df$rep)
df$mid = (df$bp.end + df$bp.start)/2
df$bp.len = (df$bp.end - df$bp.start)
df$Pi.1 = df$nuc.diversity.within_1 / df$bp.len
df$Pi.2 = df$nuc.diversity.within_2 / df$bp.len
df$dxy.1.2 = df[['dxy']] / df$bp.len

tdf = melt(df[,c('rep','mid','Pi.1','Pi.2','dxy.1.2')],id=c("rep","mid"))
pdf(file = gsub(".tmp.txt", ".plots.pdf", args[2]))
tdf$id = paste(tdf$rep, tdf$variable)
ggplot() + geom_line(data = tdf, aes(x=mid, y=value, group = id, color = variable), alpha=0.7) + geom_smooth(data = tdf, aes(x=mid, y=value, color = variable),method="loess")

ggplot() + geom_line(data = df, aes(x = mid, y = Pi.1, color = rep), alpha = 0.8) + geom_smooth(data = df, aes(x = mid, y = Pi.1,color = "green"), method="loess") + ggtitle("Population 1")

ggplot() + geom_line(data = df, aes(x = mid, y = Pi.2, color = rep), alpha = 0.8) + geom_smooth(data = df, aes(x = mid, y = Pi.2,color = "green"), method="loess") + ggtitle("Population 2")


dev.off()
      
write.table(df,args[2],row.names = FALSE, quote = FALSE, col.names = FALSE)


dat4 %>% filter(mu == 0.0000001 & recomb == 0.0000001 & N == 10000 & dom == 0.5) %>% ggplot(.) + geom_smooth(aes(x = bp.end, y = Pi.1, color = as.factor(ploidy)), linetype = "solid", method = "loess", span = 0.01, se = F, size = 0.7) + geom_smooth(aes(x = bp.end, y = Pi.2, color = as.factor(ploidy)), linetype = "dashed", method = "loess", span = 0.01, se = F, size = 0.5) + facet_grid(~s) +xlab("Position (bp)") + ylab("Diversity") + theme_bw() + scale_color_manual(name="Ploidy", values = c('red','blue'))

