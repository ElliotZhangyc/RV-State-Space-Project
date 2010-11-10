# We want to check that our inferential machine is working.  I assume
# that we have already ran <GenerateData.R> and that T was set to
# 10000.  We now look at our inferences for various portions of the
# time series, such as 1:100, 1:1000, and 1:10000.

# Select a subset of the data.
#T = 1000;
#y.data = y.true[1:T];
#z.data = z.true[1:(T+1)];

# Estimate the parametrs.
#source("CalibrateData.R");

# Output.  Comment out what you do not want.

## TO VIEW TIME SERIES ##
#jpeg("DLM4_Test_TS_1000.jpg", width=1024, height=512);
#plot(1:T, y.data, type="l");
#lines(0:T, z.data, col="grey");
#dev.off();

## FOR SUMMARY STATISTICS ##
# QoI = quantity of interest.
#range = 3000:num.samples;
#QoI = W[range];
#m.QoI = mean(QoI);
#var.QoI = var(QoI);
#print(c(m.QoI, var.QoI));

## FOR HISTOGRAMS ##
#jpeg("DLM4_mu_hist_1000.jpg", width=512, height=512);
#hist(W[range], breaks=40, prob=TRUE);
#dev.off();

## FOR AUTOCORRELATIONS ##
#the.acf = acf(phi);
#jpeg("DLM4_phi_acf_1000.jpg");
#plot(the.acf, type="l");
#dev.off();

