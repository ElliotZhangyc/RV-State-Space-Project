# We want to check that our inferential machine is working.  I assume
# that we have already ran <GenerateData.R> and that T was set to
# 10000.  We now look at our inferences for various portions of the
# time series, such as 1:100, 1:1000, and 1:10000.

# Select a subset of the data.
T = 1000;
#y.data = y.true[1:T];
#z.data = z.true[1:(T+1)];

# Estimate the parametrs.
source("CalibrateData.R");

# Output.  Comment out what you do not want.

## TO VIEW TIME SERIES ##
#jpeg("DLM3_Test_TS_1000.jpg", width=1024, height=512);
#plot(1:T, y.data, type="l");
#lines(0:T, z.data, col="grey");
#dev.off();

## FOR SUMMARY STATISTICS ##
# QoI = quantity of interest.
range = 1000:num.samples;
QoI = W[range];
m.QoI = mean(QoI);
var.QoI = var(QoI);
print(c(m.QoI, var.QoI));

#par(mfrow=c(2,3))
#par(mai=c(0.25, 0.25, 0.25, 0.25))

## FOR HISTOGRAMS ##
#jpeg("DLM3_W_hist_Test2_10000_steps_10000_MCMC.jpg", width=512, height=512);
#hist(mu[range], breaks=40, prob=TRUE, main=expression(paste("Histogram of ",mu,"|y")));
#hist(W[range], breaks=40, prob=TRUE, main=expression(paste("Histogram of ",W,"|y")));
#hist(phi[range], breaks=40, prob=TRUE, main=expression(paste("Histogram of ",phi,"|y")));
#dev.off();

#par(mai=c(0.40, 0.25, 0.60, 0.25))
## FOR AUTOCORRELATIONS ##
#the.acf.mu = acf(mu[1:num.samples], plot=TRUE, main=expression(paste("Autocorrelation of ",mu,"|y")));
#the.acf.W = acf(W[1:num.samples], plot=TRUE, main=expression(paste("Autocorrelation of ",W,"|y")));
#the.acf.phi = acf(phi[1:num.samples], plot=TRUE, main=expression(paste("Autocorrelation of ",phi,"|y")));
#jpeg("DLM3_W_acf_Test2_10000_steps_10000_MCMC.jpg");
#plot(the.acf, type="l");
#dev.off();

## FOR Z ##
#NN = 100;
NN = T;
z.post.mean = 1:(NN+1);
z.post.sd = 1:(NN+1);
for(i in 1:(NN+1)){
  z.post.mean[i] = mean(z[i,range]);
  z.post.sd[i] = sd(z[i,range]);
}

# Plot y.data along with its estimate z.post.mean.
#####jpeg("DLM4_z_first100_1000_steps_10000_MCMC.jpg", width=1024, height=512)
#plot(1:T, y.data[1:T]+0.65, type="l", col="grey");
#lines(0:T, z.post.mean);
#lines(0:T, z.post.mean+z.post.sd, col="grey");
#lines(0:T, z.post.mean-z.post.sd, col="grey");
#####dev.off();

# Plot z.data along with its estimate z.post.mean.
#jpeg("DLM4_z_first100_1000_steps_10000_MCMC.jpg", width=1024, height=512)
#plot(1:NN, y.data[1:NN], type="l", col="grey");
#lines(0:NN, z.data[1:(NN+1)], type="l", col="black");
#lines(0:NN, z.post.mean, col="blue");
#lines(0:NN, z.post.mean+z.post.sd, col="green");
#lines(0:NN, z.post.mean-z.post.sd, col="green");
#dev.off();

plot(1:T, y.data, type="l", col="grey");
lines(0:T, z.post.mean, col="black");
lines(0:T, z.post.mean+z.post.sd, col="green");
lines(0:T, z.post.mean-z.post.sd, col="green");
