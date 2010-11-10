# We want to check that our inferential machine is working.  I assume
# that we have already ran <GenerateData.R> and that T was set to
# 10000.  We now look at our inferences for various portions of the
# time series, such as 1:100, 1:1000, and 1:10000.

# Select a subset of the data.
#T = 100;
#y.data = y.true[1:T];
#z.data = z.true[1:(T+1)];
#q.data = q.true[1:T];

# Estimate the parametrs.
# source("CalibrateData.R");

# Output.  Comment out what you do not want.

## TO VIEW TIME SERIES ##
#jpeg("DLM4_Test_TS_100.jpg", width=1024, height=512);
#plot(1:T, y.data, type="l");
#lines(0:T, z.data, col="grey");
#dev.off();

## FOR SUMMARY STATISTICS ##
# QoI = quantity of interest.
range = 1000:num.samples;
QoI = mu[range];
m.QoI = mean(QoI);
var.QoI = var(QoI);
print(c(m.QoI, var.QoI));

## FOR HISTOGRAMS ##
jpeg("DLM4_mu_hist_SP500_2480_steps_3000_MCMC.jpg", width=512, height=512);
hist(mu[range], breaks=40, prob=TRUE);
dev.off();

## FOR AUTOCORRELATIONS ##
the.acf = acf(mu);
jpeg("DLM4_mu_acf_SP500_2480_steps_3000_MCMC.jpg", width=512, height=512);
plot(the.acf, type="l");
dev.off();

## FOR Z ##
#z.post.mean = 1:(T+1);
#z.post.sd = 1:(T+1);
#for(i in 1:(T+1)){
#  z.post.mean[i] = mean(z[range,i]);
#  z.post.sd[i] = sd(z[range,i]);
#}

# Plot y.data along with its estimate z.post.mean.
#####jpeg("DLM4_z_first100_1000_steps_10000_MCMC.jpg", width=1024, height=512)
#plot(1:T, y.data[1:T]+0.65, type="l", col="grey");
#lines(0:T, z.post.mean);
#lines(0:T, z.post.mean+z.post.sd, col="grey");
#lines(0:T, z.post.mean-z.post.sd, col="grey");
#####dev.off();

# Plot z.data along with its estimate z.post.mean.
#jpeg("DLM4_z_first100_1000_steps_10000_MCMC.jpg", width=1024, height=512)
#plot(0:T, z.data, type="l", col="blue");
#lines(0:T, z.post.mean);
#lines(0:T, z.post.mean+z.post.sd, col="grey");
#lines(0:T, z.post.mean-z.post.sd, col="grey");
#dev.off();

# Plot the difference between the actual z and the posterior mean of z.
#jpeg("DLM4_z_diff_1000_steps_10000_MCMC.jpg", width=1024, height=512);
#plot(0:T, z.data-z.post.mean, type="h");
#barplot(z.data-z.post.mean);
#dev.off();

## FOR q ##
#num.correct = 0;
#for(i in 1:T){
#  hist.q = hist(q[,i], breaks=c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5),
#                plot=FALSE);
#  themax = max(hist.q$counts);
#  idx = hist.q$mids%*%(themax==hist.q$counts);
#  num.correct = num.correct + 1*(q.data[i]==idx);
#}
#print(num.correct);
