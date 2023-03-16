rm(list=ls())
options(scipen = 999)
library("optparse")
 
option_list = list(
    make_option(c("-p", "--painting"), type="character", default=NULL, 
        help="Path to the ancestry painting vector.", metavar="character"),
    make_option(c("-f", "--fraction"), type="character", default=NULL, 
        help="Path to the ancestry fraction vector.", metavar="character"),
    make_option(c("-o", "--outputfilepath"), type="character", default=NULL, 
        help="output file path.", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$painting) | is.null(opt$fraction) | is.null(opt$outputfilepath)){
    print_help(opt_parser)
    print("[R-Outsourced Optimizing Script] All arguments are necessaary - arguments sent are not proper.")
    quit()
}

print("[R-Outsourced Optimizing Script] Optimizing and writing population-level estimates")

referenceMeans = as.matrix(read.table(opt$painting, row.names = 1, header = T, sep="\t", stringsAsFactors = F, check.names = F))
agg = read.table(opt$fraction, row.names = 1, header = T, sep="\t", stringsAsFactors = F, check.names = F)

regularizedRSS = function(data, par) {
  
    ## Data is the observed copying rates
    ## Parameters are the estimated fractions
    ## Find the difference between predicted and observed copying rates, and regularize
    predicted = t(referenceMeans) %*% cbind(par)
    error = sum((data - predicted)**2)
    penalty = error * (sum(par > 0) - 1) * .1
    #penalty = log(length(par)) * log(max(c(0, sum(par > 0) - 1)) + 1)
    return(error + penalty)
  
}

regularizedEstimates = apply(agg, MAR = 1, FUN = function(i) optim(par = rep(1/nrow(referenceMeans), nrow(referenceMeans)), fn = regularizedRSS, lower = rep(0, 5), upper = rep(1, 5), data = i, method = 'L-BFGS-B')$par)
regularizedEstimates = apply(t(regularizedEstimates), MAR = 1, function(i) i / sum(i))
rownames(regularizedEstimates) = rownames(referenceMeans)
estimates = regularizedEstimates
estimates = t(estimates)
write.table(x=round(estimates, 6), file=opt$outputfilepath, sep="\t", quote = F)

