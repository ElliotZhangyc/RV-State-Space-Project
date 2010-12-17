
# PLOTTING SOME STUFF

## HISTOGRAMS AND CORRELOGRAM ##

# If you want to write the plot to file.
#postscript("pairwise-compare.eps", width=6,
#           height=2, horizontal = FALSE, onefile = FALSE, paper = "special");

# Set the margins.
par(mai=c(0.40, 0.25, 0.4, 0.25));

# Set up the layout.
layoutmat = matrix(c(1,2,3), 1, 3, byrow=TRUE);
layout(layoutmat);

# Histogram of gamma
hist(gamma.gibbs[range], breaks=40, prob=TRUE,
     main=expression(paste("Histogram of ",gamma,"|y")));
# Histogram of mu
hist(mu.gibbs[range], breaks=40,
     main=expression(paste("Histogram of ", mu, "|y")));
# Correlogram of gamma and mu
ccf(gamma.gibbs, mu.gibbs,
     main=expression(paste("Correlogram of ", gamma, " and ", mu)));

# Remember to turn of the postscript device.
#dev.off();

## Ploting Standardized Gibbs Samples ##

# Set up the layout.
#layoutmat = matrix(c(1,2), 1, 1, byrow=TRUE);
#layout(layoutmat);

# Standardize the time series.
#std.gamma.gibbs = (gamma.gibbs-mean(gamma.gibbs))/sd(gamma.gibbs);
#std.mu.gibbs = (mu.gibbs-mean(mu.gibbs))/sd(mu.gibbs);

# Plot a subset
#every10 = seq(1,T,10);

# Plot the Gibbs Samples.
#plot(every10, std.gamma.gibbs[every10], type="l");
#lines(every10, std.mu.gibbs[every10+1], col="blue");
