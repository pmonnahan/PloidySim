library(rio)
library(rehh)
library(readr)
library(plyr)
library(dplyr)
library(magrittr)
library(inline)
library(parallel)
library(stringr)
library(doParallel)

# Read in arguments
args = commandArgs(trailingOnly=TRUE)

ms_outdir_path = args[1]
outfile = args[2]
files <- list.files(path = ms_outdir_path, pattern = "*.out", full.names=TRUE, recursive=FALSE)

num_cores=4
# cl = makeCluster(num_cores - 1, outfile=paste(ms_outdir_path,".log"))
cl = makeCluster(num_cores - 1, outfile="")
registerDoParallel(cl)

make_flag <- function(flag, vals, strict=FALSE) {
  # make a unix CL flag by appending one or two dashes. underscore turned into -
  dash <- "-"
  if (nchar(flag) > 1 && flag != 'seed')
    dash <- "--"
  # if a logical value is passed, e.g. T=TRUE, turn that into -T
  if (is.logical(vals)) {
    if (vals)
      return(sprintf("%s%s", dash, flag))
    # else, flag is off, don't include 
  }
  # if a underscore is present (e.g. for mspms's -random_seed) convert to dash
  flag <- gsub('_', '-', flag, fixed=TRUE) 
  sprintf("%s%s %s", dash, flag, paste(format(vals, scientific=FALSE), collapse=" "))
}

strict_make_flag <- function(flag, vals) {
  flag <- gsub('.', '-', flag, fixed=TRUE)
  sprintf("%s %s", flag, paste(format(vals, scientific=FALSE), collapse=" "))
}

make_args <- function(..., strict=FALSE) {
  args <- list(...)
  if (!strict)
    return(paste(Map(make_flag, names(args), args), collapse=" "))
  else
    return(paste(Map(strict_make_flag, names(args), args), collapse=" "))
}

#' Call MS from R
#'
#' Call Hudon's MS from R, returning a string of results (that can be parsed
#' with parse_ms()). Other coalescent simulators are supported as long as they
#' output results in MS format; specify with the argument \code{ms='mspms'}
#' (in this example, running msprime here). No argument checking occurs in this 
#' function; arguments passed to \code{...} are converted to command line 
#' arguments, e.g. \code{r=c(4, 1000)} is converted to \code{-r 4 1000} and 
#' \code{mutation.rate=1e-5} converted to \code{--mutation-rate 1e-5}.
#'
#' @param nsam number of samples (gametes) to draw
#' @param howmany how many replicates to run
#' @param cmd the command to pass to MS 
#' @param ... command line arguments as function arguments (if not using \code{cmd})
#' @param strict don't try to to guess flags by length of flag name; ".f 3" converted to -f 3, "..flag 3" converted to --flag 3
#' @param ms the ms-like executable to run, must be in \code{$PATH} or path to executable
#' 
#' @export
call_ms <- function(nsam, howmany, cmd=NULL, ..., ms="ms", strict=FALSE, verbose=TRUE) {
  func_args <- make_args(...)
  if (length(func_args) > 0 && !is.null(cmd))
    stop("specify string command line arguments through 'cmd' or arguments through '...', not both")
  if (length(func_args) > 0)
    cmd <- make_args(..., strict=strict)
  ms_cmd <- sprintf("%s %d %d %s", ms, nsam, howmany, cmd)
  if (verbose)
    message(sprintf("command: %s\n", ms_cmd))
  system(ms_cmd, intern=TRUE)
}

#' Call MS and Parse Output from R
#'
#' Call Hudon's MS from R, returning a tibble of results that has been parsed
#' by parse_ms(). This is simply a wrapper for call_ms() and parse_ms()
#'
#' @param nsam number of samples (gametes) to draw
#' @param howmany how many replicates to run
#' @param cmd the command to pass to MS 
#' @param ... command line arguments as function arguments (if not using \code{cmd})
#' @param ms the ms-like executable to run, must be in \code{$PATH} or path to executable
#
#' 
#' @export
ms <- function(nsam, howmany, cmd=NULL, ..., ms="ms", verbose=TRUE) {
  parse_ms(call_ms(nsam=nsam, howmany=howmany, cmd=cmd, ..., ms=ms, 
                   verbose=verbose))
}


#' A "tee" (in Unix sense) that writes MS's results to file
#'
#' @param x MS results @param con to save output to
#' @param pass if TRUE, pass the input directly to output for pipe operations
#' 
#' This writes a connection \code{con}, then return results (so can be used in
#' pipe).
#'
#' @export
#' @examples
#' # write MS results to 'test.ms', then pass down the pipeline
#' \dontrun{
#' res <- call_ms(10, 500, "-t 5") %>% write_tee('test.ms') %>% 
#'              parse_ms() %>% mutate(pi=map_dbl(gametes, pi))
#' }
write_tee <- function(x, con, pass=TRUE) {
  writeLines(x, con=con)
  if (pass)
    return(x)
}

#' Parse MS's key/value pairs, e.g. segsites and positions
#' returning a list of key/vals (where vals can be list too)
#' @keywords internal
parse_keyvals <- function(x) {
  keyvals <- gsub("(\\w+):\\s+(.*)", "\\1;;\\2", x, perl=TRUE)
  tmp <- strsplit(keyvals, ";;")[[1]]
  key <- tmp[1]
  vals <- as.numeric(strsplit(tmp[2], "\\s+")[[1]])
  if (length(vals) == 1)
    return(setNames(tibble::tibble(vals), key))
  else
    return(setNames(tibble::tibble(list(vals)), key))
}


#' Convert gaemtes' alleles intro matrix
#' @keywords internal
sites_matrix <- function(x) {
  do.call(rbind, lapply(x, function(y) as.integer(unlist(strsplit(y, "")))))
}

has_tree <- function(x) {
  # x are lines for a sim replicate. If there's a tree, second line
  # begins with (
  regexpr("^(\\(|\\[)", x[2]) != -1
}

has_many_trees <- function(x) regexpr("^\\[", x[2]) != -1

extract_recomb_trees <- function(x) {
  chunks <- strsplit(sub("\\[(\\d+)\\](.+)", "\\1;;;;\\2", x, perl=TRUE), ";;;;")
  lens <- as.integer(sapply(chunks, '[', 1))
  trees <- sapply(chunks, '[', 2)
  list(tibble::tibble(lens=lens, tree=trees))
}

#' Tidy a single simulation result from MS
#' @keywords internal
tidy_sim <- function(x) {
  # first element is the delimiter "//", last element is blank line (except for
  # last sim) 
  stopifnot(x[1] == "//")
  i <- 2
  keyval_i <- grep("^[a-z]+", x)
  if (has_tree(x)) {
    # for sims with recombination, we have multiple trees until segsites
    if (has_many_trees(x)) {
      # parse each tree, extracting the segment length (bp) and tree
      tree <- extract_recomb_trees(x[i:(min(keyval_i)-1)])
    } else {
      tree <- x[i]  # just one tree per locus, e.g. no recomb
    }
  }
  
  keyval_lines <- x[keyval_i]
  
  # parse keyval lines, and start creating the master tibble
  out <- dplyr::bind_cols(lapply(keyval_lines, parse_keyvals))
  
  stopifnot(nrow(out) == 1)
  if (out$segsites[1] > 0) {
    # extract gamete lines, remove empty line if there, convert gametes to matrices
    gametes <- x[(max(keyval_i)+1):length(x)]
    gametes <- gametes[nchar(gametes) > 0]
    gametes <- sites_matrix(gametes)
    out$gametes <- list(gametes)
  } else {
    out$positions <- list(numeric())
    out$gametes <- list(NULL)
  }
  
  if (has_tree(x))
    out$tree <- tree
  out
}

parse_ms <- function(x, include_seeds=FALSE) {
  cmd <- x[1]
  seeds <- as.integer(strsplit(x[2], " ")[[1]])
  x <- x[4:length(x)]  # drop first few lines of metadata
  res_grp <- cumsum(x == "//")
  res <- split(x, res_grp)
  sims_lst  <- lapply(res, tidy_sim)
  sims <- do.call(rbind, sims_lst)
  if (include_seeds)
    sims$seeds <- paste(seeds, collapse=' ')
  # sims <- bind_rows(sims_lst)  # bind_rows() fails over 1000 entries
  sims$rep <- seq_along(sims_lst)
  
  colorder <- c('rep', 'segsites', 'positions', 'gametes')
  if (include_seeds)
    colorder <- c('rep', 'seeds', 'segsites', 'positions', 'gametes')
  if ('time' %in% colnames(sims)) {
    # unpack this special column
    sims <- mutate(sims, tmrca=map_dbl(time, first),
                   ttot=map_dbl(time, nth, n=2)) 
    colorder <- c(colorder, 'tmrca', 'ttot')
  }
  if ('tree' %in% colnames(sims))
    colorder <- c(colorder, 'tree')
  out <- sims[, colorder]
  out
}

RE.input = function(ms_file, jit_amt = 1){
  bb = read_lines(ms_file)
  aa = bb %>% parse_ms()
  haps = as.numeric(str_split(bb[1]," ")[[1]][2])
  # L = as.numeric(str_split(bb[1]," ")[[1]][10])
  L=1000000
  hap_file1 = str_replace(ms_file,"[.]out",".hap1") #hap_file population 1
  hap_file2 = str_replace(ms_file,"[.]out",".hap2") #population 2 - nonsweep
  map_file = str_replace(ms_file,"[.]out",".map")
  
  write.table(aa$gametes[[1]][1:(haps/2),1:aa$segsites] + 1, hap_file1, col.names = F)
  write.table(aa$gametes[[1]][((haps/2) + 1):haps,1:aa$segsites] + 1, hap_file2, col.names = F)
  
  write.table(data.frame("snp" = seq(1,aa$segsites), "chr" = rep("1",aa$segsites), "pos" = sort(ceiling(jitter(aa$positions[[1]] * L, amount = jit_amt))), "anc" = rep(1, aa$segsites), "der" = rep(2, aa$segsites)), map_file, col.names = F, row.names = F)
  return(c(hap_file1, hap_file2, map_file))
}

calcHAP = function(HapPaths){
  invisible(capture.output(hh1 <- data2haplohh(HapPaths[1], HapPaths[3]))) #Import data
  invisible(capture.output(hh2 <- data2haplohh(HapPaths[2], HapPaths[3])))
  ss1 = scan_hh(hh1) #Calc iHH and iES
  ss2 = scan_hh(hh2)
  
  rsb = ies2rsb(ss1,ss2) #These are cross population haplotype estimates of selective sweep
  xp.ehh = ies2xpehh(ss1,ss2)
  
  colnames(rsb)[4] = "rsb.p"
  colnames(xp.ehh)[4] = "xpehh.p"
  HAP = merge(rsb, xp.ehh, by=c("CHR","POSITION"))
  return(HAP)
}

rehh_parallel = function(j, delete=TRUE) {
  x = files[j]
  params = str_split(x,"_")[[1]]
  
  paths = tryCatch({
    RE.input(x)
  },
  warning = function(w){
    paths = RE.input(x)
    return(paths)
  },
  error = function(e){
    message("Error: RE.input")
    message(x)
    return(NULL)
  })
  print(paths)
  if (!is.null(paths)){
    out = tryCatch({
      out <- calcHAP(paths)
    }, 
    warning = function(w){
      message("Warning: calcHAP")
      message(x)
      return(NULL)
    },
    error = function(e){
      message("Error: calcHAP")
      message(x)
      return(NULL)
    })} else {
      message("RE.input failure")
      message(x)
    }
  # add params
  if (!is.null(out)){
    ndat = tryCatch({
      out %>% mutate(traj_ploidy = substr(params[1], nchar(params[1]), nchar(params[1])), sim_ploidy = substr(params[2],3,3), dom = params[3], s = params[4], N = params[5], mu = params[6], recomb = params[7], sampGen = params[8], fuseGen = params[9], rep = params[10])
    }, 
    error = function(e){
      message("Error: Add params")
      })
  }
  if (delete==TRUE){
    file.remove(paths[1])
    file.remove(paths[2])
    file.remove(paths[3])
  }
  return(ndat)
}

# Call to parallel
# 
# includes <- '#include <sys/wait.h>'
# code <- 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
# wait <- cfunction(body=code, includes=includes, convention='.C')
# 
# for (i in 1:(length(files)/num_cores)){
#   j = i * num_cores
#   progress = round(i/(length(files)/num_cores), digits = 4) * 100
#   print(paste0("Running jobs: ", j - (num_cores-1), "-", j, "; (", progress, "% complete)"))
#   ndat = mclapply((j - (num_cores-1)):j, FUN = rehh_parallel, mc.cores = num_cores)
#   wait()
#   ndat = as.data.frame(do.call(rbind.fill, ndat))
#   if (i==1){
#     dat = ndat
#   } else {
#     dat = rbind.fill(ndat, dat)
#   }
# }

info = file.info(files)
empty = rownames(info[info$size == 0, ])
files = files[!files %in% empty]

dat = foreach(i = 1:(length(files)), .combine=rbind, 
        .packages=c('rio','rehh','readr','dplyr','inline','stringr','magrittr'),.export=c("calcHAP","RE.input"), .errorhandling="remove") %dopar% {
        print(paste(i," of ",length(files),"; ", round(i/length(files)*100, 2), "% complete", sep = ''))
        params = str_split(basename(files[i]),"_")[[1]]
        paths = tryCatch({
            RE.input(files[i])
          },
          error = function(e){
            message("Error: RE.input")
            message(paths)
            return(NULL)
          })
        
        if (!is.null(paths)){
          out = tryCatch({calcHAP(paths)
            },
            error = function(e){
              message("Error: calcHAP")
              message(paths[3])
              return(NULL)
            })
        }
        if (!is.null(out)){
          out %<>% mutate(traj_ploidy = params[1], sim_ploidy = params[2], dom = params[3], s = params[4], N = params[5], mu = params[6], recomb = params[7], sampGen = params[8], fuseGen = params[9], rep = params[10])
          file.remove(paths)
          out
        }  
        }
print(warnings())
stopImplicitCluster()

write.csv(dat, outfile, quote=F, row.names = F)

# FF = RE.input(ms_file)
# dat = calcHAP(FF)
# params = str_split(ms_file,"_")[[1]]
# 
# hap = do.call("rbind",lapply(files, function(x) {
#   FF = RE.input(x)
#   # apply function
#   params = str_split(x,"_")[[1]]
#   out <- calcHAP(FF)
#   out %<>% mutate(traj_ploidy = substr(params[1], nchar(params[1]), nchar(params[1])), sim_ploidy = substr(params[2],3,3), dom = params[3], s = params[4], N = params[5], mu = params[6], recomb = params[7], sampGen = params[8], fuseGen = params[9], rep = params[10])
#   return(out)
# }))