
set.seed(1)
pvals = runif(1e3, 0, 1)
pvals[1] =  NA
pvals[2] = Inf

qq_plot(pvals, full = 1) 
