# We want to check that our inferential machine is working.  

## FOR SUMMARY STATISTICS ##

# Summary statistics for mu, phi, and W.
range = mcmc$burn.in:mcmc$samples;
post = data.frame(
  m.mu = mean(mu.gibbs[range]),
  v.mu = var(mu.gibbs[range]),
  m.phi = mean(phi.gibbs[range]),
  v.phi = var(phi.gibbs[range]),
  m.W = mean(W.gibbs[range]),
  v.W = var(W.gibbs[range])
)

## FOR Z ##
z.post.mean = 1:(T+1);
z.post.sd = 1:(T+1);
for(i in 1:(T+1)){
  z.post.mean[i] = mean(z.gibbs[i,range]);
  z.post.sd[i] = sd(z.gibbs[i,range]);
}

## SUMMARY PLOTS ##

## Uncomment this to output to postscript.  Remeber to uncomment dev.off().
postscript(file=paste("plots_", trial.id, ".eps", sep=""), width=6,
           height=6, horizontal = FALSE, onefile = FALSE, paper = "special");

## SETUP THE PLOT PARAMETERS ##

# The margins.
par(mai=c(0.40, 0.25, 0.4, 0.25));
# The layout of the plot.
layoutmat = matrix(c(1,1,1,2,3,4,5,6,7), 3, 3, byrow=TRUE);
layout(layoutmat);

## CREATE PLOTS ##

## FOR Y AND Z ##
plot.length = min(T, 100);
y.range = 1:plot.length;
z.range = 1:(plot.length+1);
plot(y.range, y.data[y.range], type="l", col="black",
     main="The time series and estimates.");
lines(0:plot.length, z.post.mean[z.range], col="blue");
lines(0:plot.length, z.post.mean[z.range]+z.post.sd[z.range], col="grey");
lines(0:plot.length, z.post.mean[z.range]-z.post.sd[z.range], col="grey");
if (is.synthetic) lines(0:plot.length, z.data[z.range], col="red", lty=2);

## FOR HISTOGRAMS ##
hist(mu.gibbs[range], breaks=40, prob=TRUE,
     main=expression(paste("Histogram of ",mu,"|y")));
hist(W.gibbs[range], breaks=40, prob=TRUE,
     main=expression(paste("Histogram of ",W,"|y")));
hist(phi.gibbs[range], breaks=40, prob=TRUE,
     main=expression(paste("Histogram of ",phi,"|y")));

## FOR AUTOCORRELATIONS ##
acf(mu.gibbs, plot=TRUE,
    main=expression(paste("Autocorrelation of ",mu,"|y")));
acf(W.gibbs, plot=TRUE,
    main=expression(paste("Autocorrelation of ",W,"|y")));
acf(phi.gibbs, plot=TRUE,
    main=expression(paste("Autocorrelation of ",phi,"|y")));

## DEVICE ON/OFF ##
dev.off()

## TO MAKE VARIOUS TABLES ##

## ASIDE: I could also use write.table and then use a Bash or Perl
## script to generate my latex table.

# For the true values.
latex.table(round(true,3), paste("Tbl-true-", trial.id, sep=""),
            table.env=FALSE);

# For the seed values.
latex.table(round(seed,3), paste("Tbl-seed-", trial.id, sep=""),
            table.env=FALSE);

# For the MCMC information.
latex.table(round(mcmc,3), paste("Tbl-mcmc-", trial.id, sep=""),
            table.env=FALSE);

# For the prior parameters.
latex.table(round(prior,3), paste("Tbl-prior-", trial.id, sep=""),
            table.env=FALSE);

# For the posterior statisics.
latex.table(round(post,3), paste("Tbl-post-", trial.id, sep=""),
            table.env=FALSE);
